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
} PLOTGROUP;

double headroom = 0.24;
double scalepos = 60.0;

/* Backing pixmap for drawing area */
static GdkPixmap *pixmap = NULL;

extern void file_selector (char *msg, char *startdir, int action, 
                           gpointer data);

int ps_print_plots (const char *fname, int flag, gpointer data);
static int five_numbers (gpointer data);
static void plot_to_xpm (gpointer data);

#ifdef G_OS_WIN32
static int cb_copy_image (gpointer data);
#endif

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
                      SAVE_BOXPLOT_EPS, data);
    }
    else if (key->keyval == GDK_p) {  
	five_numbers(data);
    }
    else if (key->keyval == GDK_b) {  
	plot_to_xpm(data);
    }
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
        file_selector("Save boxplot file", paths.userdir, 
                      SAVE_BOXPLOT_EPS, ptr);
    else if (!strcmp(item, "Save as PS...")) 
        file_selector("Save boxplot file", paths.userdir, 
                      SAVE_BOXPLOT_PS, ptr);
#ifdef G_OS_WIN32
    else if (!strcmp(item, "Copy to clipboard"))
	cb_copy_image(data);
#endif
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
#endif	
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
				 "Helvetica", 12,
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
    double start = scalepos + plotgrp->width * (headroom / 2.0);
    double xrange = (1.0 - headroom) * plotgrp->width - scalepos;
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
gtk_boxplot_yscale (PLOTGROUP *grp, GtkPlotPC *pc)
{
    double points[4];
    double top = (headroom / 2.0) * grp->height;
    double bottom = (1.0 - headroom / 2.0) * grp->height;
    char numstr[16];
    GdkGC *gc = NULL;

    if (pc == NULL) gc = grp->window->style->bg_gc[GTK_STATE_SELECTED];

    /* draw vertical line */
    points[0] = points[2] = scalepos;
    points[1] = top;
    points[3] = bottom;
    draw_line (points, grp->area, gc, pc);

    /* draw backticks top and bottom */
    points[2] = points[0] - 5;
    points[1] = points[3] = top;
    draw_line (points, grp->area, gc, pc);
    points[1] = points[3] = bottom;
    draw_line (points, grp->area, gc, pc);

    /* draw backtick at middle */
    points[1] = points[3] = top + (bottom - top) / 2.0;
    draw_line (points, grp->area, gc, pc);
    
    /* mark max and min values on scale */
    sprintf(numstr, "%.4g", grp->gmax);
    setup_text (grp->area, gc, pc, numstr, scalepos - 10, top, 
		GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", grp->gmin);
    setup_text (grp->area, gc, pc, numstr, scalepos - 10, bottom, 
		GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", (grp->gmax + grp->gmin) / 2.0);
    setup_text (grp->area, gc, pc, numstr, scalepos - 10, 
		top + (bottom - top) / 2.0, GTK_JUSTIFY_RIGHT);
}

/* ............................................................. */

static void 
gtk_area_boxplot (BOXPLOT *plot, GtkWidget *area, 
		  GtkStyle *style, GtkPlotPC *pc,
		  int height, double boxwidth, 
		  double gmax, double gmin)
{
    double points[4];
    double ybase = height * headroom / 2.0;
    double xcenter = plot->xbase + boxwidth / 2.0;
    double scale = (1.0 - headroom) * height / (gmax - gmin);
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
		height * (1.0 - headroom/4.0), GTK_JUSTIFY_CENTER);
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
		       GTK_SIGNAL_FUNC(configure_event), NULL);

    gtk_signal_connect(GTK_OBJECT(grp->area), "expose_event",
		       GTK_SIGNAL_FUNC(expose_event), NULL);

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
	ps->page_width = ps->pc.width = pscale * grp->width;
	ps->page_height = ps->pc.height = pscale * grp->height;
    }

    if (!psinit (ps)) return 1;
    gtk_psfont_init();

    for (i=0; i<grp->nplots; i++)
	gtk_area_boxplot (&grp->plots[i], 
			  NULL, NULL, &ps->pc, 
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

    pprintf(&prn, "Five-number summar%s\n\n"
	    "%22s%12s%12s%12s%12s\n",
	    (grp->nplots > 1)? "ies" : "y",
	    "minimum", "Q1", "median", "Q3", "maxmimum");

    for (i=0; i<grp->nplots; i++) {
	pprintf(&prn, "%-10s%12g%12g%12g%12g%12g\n",
		grp->plots[i].varname, grp->plots[i].min, 
		grp->plots[i].lq, grp->plots[i].median,
		grp->plots[i].uq, grp->plots[i].max);
    }
    (void) view_buffer(&prn, 78, 240, "gretl: 5 numbers", BXPLOT,
                       view_items);

    return 0;
}

/* ............................................................. */

int boxplots (const int *list, double **pZ, const DATAINFO *pdinfo)
{
    int i, n = pdinfo->t2 - pdinfo->t1 + 1;
    double *x;
    PLOTGROUP *plotgrp;
    int width = 600, height = 450;

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

    plotgrp->height = height;
    plotgrp->width = width;
    place_plots (plotgrp);

    if (make_area(plotgrp) == NULL) return 1;
    plotgrp->popup = NULL;
    gtk_widget_show(plotgrp->window);

    for (i=0; i<plotgrp->nplots; i++)
	gtk_area_boxplot (&plotgrp->plots[i], 
			  plotgrp->area, plotgrp->window->style, NULL, 
			  height, plotgrp->boxwidth, 
			  plotgrp->gmax, plotgrp->gmin);
    
    gtk_boxplot_yscale(plotgrp, NULL);

    return 0;
}

static void plot_to_xpm (gpointer data)
{    
    PLOTGROUP *grp = (PLOTGROUP *) data;
    GdkImage *image;
    int i, j;
    guint32 pixel, white_pixel;
    FILE *fp;
    
    fp = fopen("boxtest.xpm", "w");
    if (fp == NULL) {
	errbox ("Couldn't open XPM file for writing");
	return;
    }

    fprintf(fp, "/* XPM */\n"
	    "static char *boxplot[] = {\n"
	    "/* width height ncolors chars_per_pixel */\n"
	    "\"%d %d 2 1\",\n"
	    "/* colors */\n"
	    "\"  c white\",\n"
	    "\". c black\",\n"
	    "/* pixels */\n", grp->width, grp->height);

    image = gdk_image_get (grp->area->window, 0, 0, 
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

}

#ifdef G_OS_WIN32

static int cb_copy_image (gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    GdkImage *image;
    int i, j;
    guint32 pixel, white_pixel;
    
    int nSizeDIB = 0;
    int nSizeLine = 0; /* DIB lines are 32 bit aligned */

    HANDLE hDIB;
    BOOL bRet;

    image = gdk_image_get(grp->area->window, 0, 0, 
			  grp->width, grp->height);

    white_pixel = pow(2, image->depth) - 1;

    /* allocate room for DIB */
    nSizeLine = ((grp->width*3-1)/4+1)*4;
    nSizeDIB = nSizeLine * grp->height + sizeof (BITMAPINFOHEADER);
  
    hDIB = GlobalAlloc (GMEM_MOVEABLE | GMEM_DDESHARE, nSizeDIB);
    if (NULL == hDIB) {
	    errbox ("Failed to allocate DIB");
	    bRet = FALSE;
    }

    /* fill header info */
    if (bRet) {
	BITMAPINFOHEADER *pInfo;
      
	bRet = FALSE;
	pInfo = GlobalLock (hDIB);
	if (pInfo) {
	    pInfo->biSize   = sizeof(BITMAPINFOHEADER);
	    pInfo->biWidth  = grp->width;
	    pInfo->biHeight = -grp->height; /* top-down */
	    pInfo->biPlanes = 1;
	    pInfo->biBitCount = 1;
	    pInfo->biCompression = BI_RGB; /* none */
	    pInfo->biSizeImage = 0; /* not calculated/needed */
	    pInfo->biXPelsPerMeter =
		pInfo->biYPelsPerMeter = 0;
	    /* color map size */
	    pInfo->biClrUsed = 0;
	    pInfo->biClrImportant = 0; /* all */
          
	    GlobalUnlock (hDIB);
	    bRet = TRUE;
	} else
	    errbox("Failed to lock DIB Header");
    }
  
    /* copy data to DIB */
    if (bRet) {
	unsigned char *pData;
      
	bRet = FALSE;
	pData = GlobalLock (hDIB);
      
	if (pData) {

	    /* calculate real offset */
	    pData += sizeof(BITMAPINFOHEADER);

	    for (i=0; i<grp->height; i++) {
		for (j=0; j<grp->width; j++) {
		    pixel = gdk_image_get_pixel(image, j, i);
		    if (pixel == white_pixel) *pData++ = 0;
		    else *pData++ = 1;
		}
	    }	    
          
	    bRet = TRUE;
          
	    GlobalUnlock (hDIB);
	} /* (pData) */
	else
	    errbox("Failed to lock DIB Data");
    } /* copy data to DIB */
  
    /* copy DIB to ClipBoard */
    if (bRet) {      
	if (!OpenClipboard (NULL)) {
	    errbox ("Cannot open the Clipboard!");
	    bRet = FALSE;
	} else {
	    if (bRet && !EmptyClipboard ()) {
		errbox ("Cannot empty the Clipboard");
		bRet = FALSE;
	    }
	    if (bRet) {
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
  
    return bRet;
} 

#endif /* G_OS_WIN32 */
