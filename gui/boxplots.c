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
#endif

typedef struct {
    double median;
    double uq, lq;
    double max, min;
    double xbase;
    char varname[9];
} BOXPLOT;

typedef struct {
    int nplots;
    BOXPLOT *plots;
    int winwidth, winheight;
    double boxwidth;
    double gmax, gmin;
} PLOTGROUP;

double headroom = 0.24;
double scalepos = 60.0;

/* Backing pixmap for drawing area */
static GdkPixmap *pixmap = NULL;

extern void file_selector (char *msg, char *startdir, int action, 
                           gpointer data);

int ps_print_plots (const char *fname, gpointer data);

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
		    widget->style->fg_gc[GTK_WIDGET_STATE(widget)],
		    pixmap,
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
        file_selector("Save boxplot file", paths.userdir, 
                      SAVE_BOXPLOT, data);
    }
    else if (key->keyval == GDK_p) {  /* temporary, for testing */
	ps_print_plots("foo.eps", data);
    }
    return TRUE;
}

/* ............................................................. */

static void
setup_text (GtkWidget *area, GdkGC *gc, GtkPlotPC *pc, char *text, 
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
				 "Helvetica", 11,
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
}

/* ............................................................. */

static void 
draw_line (double *points, GtkWidget *area, GdkGC *gc, GtkPlotPC *pc)
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
    double start = scalepos + plotgrp->winwidth * (headroom / 2.0);
    double xrange = (1.0 - headroom) * plotgrp->winwidth - scalepos;
    double boxwidth = xrange / (2.0 * plotgrp->nplots - 1.0);

    plotgrp->boxwidth = boxwidth;

    /* divide horizontal space between plots; also get global
       maximum and minimum y-values */

    plotgrp->gmax = plotgrp->plots[0].max;
    plotgrp->gmin = plotgrp->plots[0].min;

    for (i=0; i<plotgrp->nplots; i++) {
	plotgrp->plots[i].xbase = start + (2.0 * (i+1.0) - 2.0) * boxwidth;
	if (i > 0) {
	    if (plotgrp->plots[i].max > plotgrp->gmax) 
		plotgrp->gmax = plotgrp->plots[i].max;
	    if (plotgrp->plots[i].min < plotgrp->gmin) 
		plotgrp->gmin = plotgrp->plots[i].min;
	}	
    }
}

/* ............................................................. */

static void 
gtk_boxplot_yscale (GtkWidget *area, GtkStyle *style, GtkPlotPC *pc,
		    int winheight, double gmax, double gmin)
{
    double points[4];
    double top = (headroom / 2.0) * winheight;
    double bottom = (1.0 - headroom / 2.0) * winheight;
    char numstr[16];
    GdkGC *gc = NULL;

    if (pc == NULL) gc = style->bg_gc[GTK_STATE_SELECTED];

    /* draw vertical line */
    points[0] = points[2] = scalepos;
    points[1] = top;
    points[3] = bottom;
    draw_line (points, area, gc, pc);

    /* draw backticks top and bottom */
    points[2] = points[0] - 5;
    points[1] = points[3] = top;
    draw_line (points, area, gc, pc);
    points[1] = points[3] = bottom;
    draw_line (points, area, gc, pc);

    /* draw backtick at middle */
    points[1] = points[3] = top + (bottom - top) / 2.0;
    draw_line (points, area, gc, pc);
    
    /* mark max and min values on scale */
    sprintf(numstr, "%.4g", gmax);
    setup_text (area, gc, pc, numstr, scalepos - 10, top, GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", gmin);
    setup_text (area, gc, pc, numstr, scalepos - 10, bottom, GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", (gmax + gmin) / 2.0);
    setup_text (area, gc, pc, numstr, scalepos - 10, 
		top + (bottom - top) / 2.0, GTK_JUSTIFY_RIGHT);
}

/* ............................................................. */

static void 
gtk_area_boxplot (BOXPLOT *plot, GtkWidget *area, 
		  GtkStyle *style, GtkPlotPC *pc,
		  int winheight, double boxwidth, 
		  double gmax, double gmin)
{
    double points[4];
    double ybase = winheight * headroom / 2.0;
    double xcenter = plot->xbase + boxwidth / 2.0;
    double scale = (1.0 - headroom) * winheight / (gmax - gmin);
    double median, uq, lq, maxval, minval;
    GdkRectangle rect;
    GdkGC *gc = NULL;
    GdkGC *whitegc = NULL;

    if (pc == NULL) {
	gc = style->bg_gc[GTK_STATE_SELECTED];
	whitegc = style->fg_gc[GTK_STATE_SELECTED];
    }

    median = ybase + (gmax - plot->median) * scale;
    uq = ybase + (gmax - plot->uq) * scale;
    lq = ybase + (gmax - plot->lq) * scale;
    maxval = ybase + (gmax - plot->max) * scale;
    minval = ybase + (gmax - plot->min) * scale;

    /* draw inter-quartile box */
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
			    TRUE, /* filled ? */
			    plot->xbase, 
			    uq, 
			    boxwidth, 
			    lq - uq);
	gtk_widget_draw (area, &rect);
    }

    /* draw line at median */
    points[0] = plot->xbase;
    points[1] = points[3] = median;
    points[2] = plot->xbase + boxwidth;
    draw_line(points, area, whitegc, pc);

    /* draw line to maximum value */
    points[0] = points[2] = xcenter;
    points[1] = maxval;
    points[3] = uq;
    draw_line(points, area, gc, pc);

    /* plus a little crossbar */
    points[0] = xcenter - 5.0;
    points[2] = xcenter + 5.0;
    points[1] = points[3] = maxval;
    draw_line(points, area, gc, pc);

    /* draw line to minimum value */
    points[0] = points[2] = xcenter;
    points[1] = lq;
    points[3] = minval;
    draw_line(points, area, gc, pc);

    /* plus a little crossbar */
    points[0] = xcenter - 5.0;
    points[2] = xcenter + 5.0;
    points[1] = points[3] = minval;
    draw_line(points, area, gc, pc);

    /* write name of variable beneath */
    setup_text (area, gc, pc, plot->varname, xcenter, 
		winheight * (1.0 - headroom/4.0), GTK_JUSTIFY_CENTER);
}

/* ............................................................. */

static void 
destroy_boxplots (GtkWidget *w, gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;

    free(grp->plots);
    free(grp);
}

/* ............................................................. */

static GtkWidget *
make_area (GtkWidget **area, PLOTGROUP *plotgrp, int winwidth, int winheight)
{
    GtkWidget *boxwin;

    boxwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(boxwin), "gretl: boxplots");

    /* Create the drawing area */
    *area = gtk_drawing_area_new ();

    gtk_signal_connect(GTK_OBJECT(*area), "configure_event",
		       GTK_SIGNAL_FUNC(configure_event), NULL);

    gtk_signal_connect(GTK_OBJECT(*area), "expose_event",
		       GTK_SIGNAL_FUNC(expose_event), NULL);

    gtk_signal_connect(GTK_OBJECT(boxwin), "key_press_event", 
		       GTK_SIGNAL_FUNC(box_key_handler), plotgrp);

    gtk_signal_connect(GTK_OBJECT(boxwin), "destroy",
		       GTK_SIGNAL_FUNC (destroy_boxplots), plotgrp);

    gtk_drawing_area_size (GTK_DRAWING_AREA(*area), winwidth, winheight); 
    gtk_widget_show (*area);

    gtk_container_add(GTK_CONTAINER(boxwin), *area);

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

/* At this point we need to borrow a few functions from gtkplotps.c,
   which is part of gtkextra.  These functions are declared as
   static in that context, but we want direct access to them */

#include "plotps.c"

/* ............................................................. */

int ps_print_plots (const char *fname, gpointer data) 
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    GtkPlotPS *ps;
    int i;
    gdouble pscale = 0.8;

    ps = GTK_PLOT_PS(gtk_plot_ps_new (fname, 
				      GTK_PLOT_PORTRAIT, 
				      TRUE, /* epsflag */
				      GTK_PLOT_LETTER, 
				      pscale, pscale)
		     );

    if (ps == NULL) return 1;

    ps->page_width = ps->pc.width = pscale * grp->winwidth;
    ps->page_height = ps->pc.height = pscale * grp->winheight;

    if (!psinit (ps)) return 1;
    gtk_psfont_init();

    for (i=0; i<grp->nplots; i++)
	gtk_area_boxplot (&grp->plots[i], 
			  NULL, NULL, &ps->pc, 
			  grp->winheight, grp->boxwidth, 
			  grp->gmax, grp->gmin);
    
    gtk_boxplot_yscale (NULL, NULL, &ps->pc, 
			grp->winheight,
			grp->gmax, grp->gmin); 

    psleave (ps);
    gtk_object_destroy(GTK_OBJECT(ps));
    gtk_psfont_unref();

    return 0;
}

/* ............................................................. */

int boxplots (const int *list, double **pZ, const DATAINFO *pdinfo)
{
    int i, n = pdinfo->t2 - pdinfo->t1 + 1;
    double *x;
    PLOTGROUP *plotgrp;
    GtkWidget *boxwin, *area;
    int winwidth = 600, winheight = 450;

    x = mymalloc(n * sizeof *x);
    if (x == NULL) return 1;

    plotgrp = mymalloc(sizeof *plotgrp);
    if (plotgrp == NULL) return 1;

    plotgrp->nplots = list[0];
    plotgrp->plots = mymalloc (plotgrp->nplots * sizeof(BOXPLOT));
    if (plotgrp->plots == NULL) {
	free(plotgrp);
	return 1;
    }

    for (i=0; i<plotgrp->nplots; i++) {
	n = ztox(list[i+1], x, pdinfo, *pZ);
	if (n < 2) {
	    fprintf(stderr, "Inadequate sample range.\n");
	    free(x);
	    free(plotgrp->plots);
	    return 1;
	}
	qsort(x, n, sizeof *x, compare_doubles);
	plotgrp->plots[i].min = x[0];
	plotgrp->plots[i].max = x[n-1];
	quartiles(x, n, &plotgrp->plots[i]);
	strcpy(plotgrp->plots[i].varname, pdinfo->varname[list[i+1]]);
    }
    free(x);

    plotgrp->winheight = winheight;
    plotgrp->winwidth = winwidth;
    place_plots (plotgrp);

    boxwin = make_area(&area, plotgrp, winwidth, winheight);
    gtk_widget_show(boxwin);

    for (i=0; i<plotgrp->nplots; i++)
	gtk_area_boxplot (&plotgrp->plots[i], 
			  area, boxwin->style, NULL, 
			  winheight, plotgrp->boxwidth, 
			  plotgrp->gmax, plotgrp->gmin);
    
    gtk_boxplot_yscale(area, boxwin->style, NULL, 
		       winheight, plotgrp->gmax, plotgrp->gmin);

    return 0;
}




