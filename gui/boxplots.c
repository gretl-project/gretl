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
/* #include <gdk-pixbuf/gdk-pixbuf.h> */

typedef struct {
    double median;
    double uq, lq;
    double max, min;
    double xbase, boxwidth;
    char varname[9];
} BOXPLOT;

double headroom = 0.24;
double scalepos = 60.0;

/* Backing pixmap for drawing area */
static GdkPixmap *pixmap = NULL;

/* ............................................................. */

/* Create a new backing pixmap of the appropriate size */
static gint
configure_event (GtkWidget *widget, GdkEventConfigure *event)
{
    if (pixmap)
	gdk_pixmap_unref(pixmap);

    pixmap = gdk_pixmap_new(widget->window,
			    widget->allocation.width,
			    widget->allocation.height,
			    -1);

    gdk_draw_rectangle (pixmap,
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
expose_event (GtkWidget *widget, GdkEventExpose *event)
{
    gdk_draw_pixmap(widget->window,
		    widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		    pixmap,
		    event->area.x, event->area.y,
		    event->area.x, event->area.y,
		    event->area.width, event->area.height);

    return FALSE;
}

/* ............................................................. */

static void
setup_text (GtkWidget *area, GdkGC *gc, char *text, 
	    double x, double y, unsigned just)
{
    GdkRectangle rect; 
    size_t len = strlen(text);
    int cw = gdk_char_width(fixed_font, 'X');

    if (just == GTK_ANCHOR_N) {
	x -= len * cw / 2.0;
    }
    else if (just == GTK_ANCHOR_E) {
	int ch = gdk_char_height(fixed_font, '1');

	y += ch / 2.0;
	x -= len * cw;
    }

    rect.x = x;
    rect.y = y;
    rect.width = 80;
    rect.height = 10;
    
    gdk_draw_string (pixmap, fixed_font, gc, x, y, text);
    gtk_widget_draw (area, &rect);
}

/* ............................................................. */

static void 
draw_line (double *points, GtkWidget *area, GdkGC *gc)
{
    GdkRectangle rect; 

    rect.x = points[0];
    rect.y = points[1];
    rect.width = points[2] - points[1] + 1;
    rect.height = points[3] - points[1] + 1;
  
    gdk_draw_line (pixmap, gc, points[0], points[1], points[2], points[3]);
    gtk_widget_draw (area, &rect);
}

/* ............................................................. */

static void 
place_plots (BOXPLOT **plotgrp, double *gmax, double *gmin,
	     int winwidth)
     /* calculate placement of plots based on the number of them
	and their respective max and min values */
{
    int i, n = 0;
    double start = scalepos + winwidth * (headroom / 2.0);
    double xrange = (1.0 - headroom) * winwidth - scalepos;
    double boxwidth;

    /* determine the number of boxplots */
    while (plotgrp[n] != NULL) n++;

    boxwidth = xrange / (2.0 * n - 1.0);

    /* divide horizontal space between plots; also get global
       maximum and minimum y-values */

    *gmax = (plotgrp[0])->max;
    *gmin = (plotgrp[0])->min;

    for (i=0; i<n; i++) {
	(plotgrp[i])->boxwidth = boxwidth;
	(plotgrp[i])->xbase = start + (2.0 * (i+1.0) - 2.0) * boxwidth;
	if (i > 0) {
	    if ((plotgrp[i])->max > *gmax) *gmax = (plotgrp[i])->max;
	    if ((plotgrp[i])->min < *gmin) *gmin = (plotgrp[i])->min;
	}	
    }
}

/* ............................................................. */

static void 
free_plotgrp (BOXPLOT **plotgrp)
{
    int i, n = 0;

    while (plotgrp[n] != NULL) n++;
    for (i=0; i<n; i++) free(plotgrp[i]);
    free(plotgrp);
}

/* ............................................................. */

static void 
gtk_boxplot_yscale (GtkWidget *area, GtkStyle *style, 
		    double gmax, double gmin, int winheight)
{
    double points[4];
    double top = (headroom / 2.0) * winheight;
    double bottom = (1.0 - headroom / 2.0) * winheight;
    char numstr[16];
    GdkGC *gc = style->bg_gc[GTK_STATE_SELECTED];

    /* draw vertical line */
    points[0] = points[2] = scalepos;
    points[1] = top;
    points[3] = bottom;
    draw_line (points, area, gc);

    /* draw backticks top and bottom */
    points[2] = points[0] - 5;
    points[1] = points[3] = top;
    draw_line (points, area, gc);
    points[1] = points[3] = bottom;
    draw_line (points, area, gc);

    /* draw backtick at middle */
    points[1] = points[3] = top + (bottom - top) / 2.0;
    draw_line (points, area, gc);
    
    /* mark max and min values on scale */
    sprintf(numstr, "%.4g", gmax);
    setup_text(area, gc, numstr, scalepos - 10, top, GTK_ANCHOR_E);
    sprintf(numstr, "%.4g", gmin);
    setup_text (area, gc, numstr, scalepos - 10, bottom, GTK_ANCHOR_E);
    sprintf(numstr, "%.4g", (gmax + gmin) / 2.0);
    setup_text (area, gc, numstr, scalepos - 10, 
		top + (bottom - top) / 2.0, GTK_ANCHOR_E);
}

/* ............................................................. */

static void 
gtk_area_boxplot (BOXPLOT *plot, double gmax, double gmin,
		  GtkWidget *area, GtkStyle *style, 
		  int winheight)
{
    double points[4];
    double ybase = winheight * headroom / 2.0;
    double xcenter = plot->xbase + plot->boxwidth / 2.0;
    double scale = (1.0 - headroom) * winheight / (gmax - gmin);
    double median, uq, lq, maxval, minval;
    GdkRectangle rect;
    GdkGC *gc = style->bg_gc[GTK_STATE_SELECTED];
    GdkGC *whitegc = style->fg_gc[GTK_STATE_SELECTED];

    median = ybase + (gmax - plot->median) * scale;
    uq = ybase + (gmax - plot->uq) * scale;
    lq = ybase + (gmax - plot->lq) * scale;
    maxval = ybase + (gmax - plot->max) * scale;
    minval = ybase + (gmax - plot->min) * scale;

    /* draw inter-quartile box */
    rect.x = plot->xbase;
    rect.y = uq;
    rect.width = plot->boxwidth;
    rect.height = uq - lq; 
    gdk_draw_rectangle (pixmap, 
			gc, 
			TRUE, /* filled ? */
			plot->xbase, 
			uq, 
			plot->boxwidth, 
			lq - uq);
    gtk_widget_draw (area, &rect);

    /* draw line at median */
    points[0] = plot->xbase;
    points[1] = points[3] = median;
    points[2] = plot->xbase + plot->boxwidth;
    draw_line(points, area, whitegc);

    /* draw line to maximum value */
    points[0] = points[2] = xcenter;
    points[1] = maxval;
    points[3] = uq;
    draw_line(points, area, gc);

    /* plus a little crossbar */
    points[0] = xcenter - 5.0;
    points[2] = xcenter + 5.0;
    points[1] = points[3] = maxval;
    draw_line(points, area, gc);

    /* draw line to minimum value */
    points[0] = points[2] = xcenter;
    points[1] = lq;
    points[3] = minval;
    draw_line(points, area, gc);

    /* plus a little crossbar */
    points[0] = xcenter - 5.0;
    points[2] = xcenter + 5.0;
    points[1] = points[3] = minval;
    draw_line(points, area, gc);

    /* write name of variable beneath */
    setup_text (area, gc, plot->varname, xcenter, 
		winheight * (1.0 - headroom/4.0), GTK_ANCHOR_N);
}

/* ............................................................. */

static GtkWidget *
make_area (GtkWidget **area, int winwidth, int winheight)
{
    GtkWidget *boxwin;
    GtkWidget *vbox;
    GtkWidget *frame;

    boxwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(boxwin), "gretl: boxplots");

    vbox = gtk_vbox_new (FALSE, 4);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 4);
    gtk_widget_show (vbox);

    /* Create the drawing area */
    *area = gtk_drawing_area_new ();

    gtk_signal_connect(GTK_OBJECT(*area), "configure_event",
		       GTK_SIGNAL_FUNC(configure_event), NULL);

    gtk_signal_connect(GTK_OBJECT(*area), "expose_event",
		       GTK_SIGNAL_FUNC(expose_event), NULL);

    gtk_drawing_area_size (GTK_DRAWING_AREA(*area), winwidth, winheight); 
    gtk_widget_show (*area);

    frame = gtk_frame_new (NULL);
    gtk_frame_set_shadow_type (GTK_FRAME (frame), GTK_SHADOW_IN);
    gtk_box_pack_start (GTK_BOX (vbox), frame, TRUE, TRUE, 0);
    gtk_widget_show (frame);

    gtk_container_add (GTK_CONTAINER (frame), *area);

    gtk_container_add(GTK_CONTAINER(boxwin), vbox);

    return boxwin;
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

static int 
compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da > *db) - (*da < *db);
}

/* ............................................................. */

int boxplots (const int *list, double **pZ, const DATAINFO *pdinfo)
{
    int i, n = pdinfo->t2 - pdinfo->t1 + 1;
    int nplots = list[0];
    double gmax, gmin, *x;
    BOXPLOT **plotgrp;
    GtkWidget *boxwin, *area;
    int winwidth = 600, winheight = 450;

    x = mymalloc(n * sizeof *x);
    if (x == NULL) return 1;

    plotgrp = mymalloc ((nplots + 1) * sizeof *plotgrp);
    if (plotgrp == NULL) return 1;
    for (i=0; i<nplots; i++) {
	plotgrp[i] = mymalloc (sizeof **plotgrp);
	if (plotgrp[i] == NULL) return 1;
    }
    plotgrp[nplots] = NULL; 

    for (i=0; i<nplots; i++) {
	n = ztox(list[i+1], x, pdinfo, *pZ);
	if (n < 2) {
	    fprintf(stderr, "Inadequate sample range.\n");
	    free(x);
	    free_plotgrp(plotgrp);
	    return 1;
	}
	qsort(x, n, sizeof *x, compare_doubles);
	(plotgrp[i])->min = x[0];
	(plotgrp[i])->max = x[n-1];
	quartiles(x, n, plotgrp[i]);
	strcpy((plotgrp[i])->varname, pdinfo->varname[list[i+1]]);
    }
    free(x);

    place_plots (plotgrp, &gmax, &gmin, winwidth);

    boxwin = make_area(&area, winwidth, winheight);
    gtk_widget_show(boxwin);

    for (i=0; i<nplots; i++)
	gtk_area_boxplot (plotgrp[i], gmax, gmin, area, 
			  boxwin->style, winheight);
    
    gtk_boxplot_yscale(area, boxwin->style, gmax, gmin, winheight);

    free_plotgrp(plotgrp);

    return 0;
}

#ifdef notdef
void grab_buffer (GtkWidget *area)
{
    GdkPixbuf *pbuf;
    GdkColormap *cmap;

    cmap = gdk_drawable_get_colormap (GTK_DRAWING_AREA(area));

    pbuf = gdk_pixbuf_get_from_drawable (NULL,
					 pixmap,
					 cmap,
					 0, 0,
					 0, 0,
					 -1, -1);

    fprintf(stderr, "pbuf is at %p\n", (void *) pbuf);
}
#endif
