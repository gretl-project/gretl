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
#include <gtkextra/gtkextra.h>

GdkColor black, white;

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
setup_text (GtkWidget *area, GdkGC *gc, GtkPlotPC *pc, char *text, 
	    double x, double y, unsigned just)
{
    GdkRectangle rect; 
    size_t len;
    int cw;

    if (pc != NULL) {
	gtk_plot_pc_draw_string (pc,
				 x, y + 4, /* alignment bodge */
				 0,
				 &black, &white,
				 FALSE,
				 0, 0, 0, 0,
				 "Courier", 11,
				 just,
				 text);
    } else {
	len = strlen(text);
	cw = gdk_char_width(fixed_font, 'X');

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
    GdkRectangle rect; 

    if (pc != NULL) {
	gtk_plot_pc_draw_line (pc, points[0], points[1], points[2], points[3]); 
    } else {
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
gtk_boxplot_yscale (GtkWidget *area, GtkStyle *style, GtkPlotPC *pc,
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
    setup_text(area, gc, pc, numstr, scalepos - 10, top, GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", gmin);
    setup_text (area, gc, pc, numstr, scalepos - 10, bottom, GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", (gmax + gmin) / 2.0);
    setup_text (area, gc, pc, numstr, scalepos - 10, 
		top + (bottom - top) / 2.0, GTK_JUSTIFY_RIGHT);
}

/* ............................................................. */

static void 
gtk_area_boxplot (BOXPLOT *plot, double gmax, double gmin,
		  GtkWidget *area, GtkStyle *style, GtkPlotPC *pc,
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

    if (pc != NULL) 
	gtk_plot_pc_draw_rectangle (pc,
				    FALSE,
				    plot->xbase, uq,
				    plot->boxwidth,
				    lq - uq);
    else {
	gdk_draw_rectangle (pixmap, 
			    gc, 
			    TRUE, /* filled ? */
			    plot->xbase, 
			    uq, 
			    plot->boxwidth, 
			    lq - uq);
	gtk_widget_draw (area, &rect);
    }

    /* draw line at median */
    points[0] = plot->xbase;
    points[1] = points[3] = median;
    points[2] = plot->xbase + plot->boxwidth;
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

/* functions borrowed from gtkplotps.c in gtkextra.  They are declared
   static there, but we want direct access to them. */

static void ps_reencode_font (FILE *file, char *fontname);

static void psleave (GtkPlotPS *ps)
{
    fprintf(ps->psfile, "showpage\n");
    fprintf(ps->psfile, "%%%%Trailer\n");
    fprintf(ps->psfile, "%%%%EOF\n");
    fclose(ps->psfile);
}

static gboolean psinit (GtkPlotPS *ps)
{
    time_t now;
    FILE *psout;

    now = time(NULL);

    psout = ps->psfile;

    if ((psout = fopen(ps->psname, "w")) == NULL){
       g_warning("ERROR: Cannot open file: %s", ps->psname); 
       return FALSE;
    }

    ps->psfile = psout;

    if (ps->epsflag)
       fprintf (psout, "%%!PS-Adobe-2.0 PCF-2.0\n");
    else
       fprintf (psout, "%%!PS-Adobe-2.0\n");

    fprintf (psout,
             "%%%%Title: %s\n"
             "%%%%Creator: %s v%s Copyright (c) 1999 Adrian E. Feiguin\n"
             "%%%%CreationDate: %s"
             "%%%%Magnification: 1.0000\n",
             ps->psname,
             "GtkPlot", "3.x",
             ctime (&now));

    if(ps->orientation == GTK_PLOT_PORTRAIT)
             fprintf(psout,"%%%%Orientation: Portrait\n");
    else
             fprintf(psout,"%%%%Orientation: Landscape\n");

    if(ps->epsflag)
          fprintf (psout,
                   "%%%%BoundingBox: 0 0 %d %d\n"
                   "%%%%Pages: 1\n"
                   "%%%%EndComments\n",
                   ps->page_width,
                   ps->page_height);


    fprintf (psout,
             "/cp {closepath} bind def\n"
             "/c {curveto} bind def\n"
             "/f {fill} bind def\n"
             "/a {arc} bind def\n"
             "/ef {eofill} bind def\n"
             "/ex {exch} bind def\n"
             "/gr {grestore} bind def\n"
             "/gs {gsave} bind def\n"
             "/sa {save} bind def\n"
             "/rs {restore} bind def\n"
             "/l {lineto} bind def\n"
             "/m {moveto} bind def\n"
             "/rm {rmoveto} bind def\n"
             "/n {newpath} bind def\n"
             "/s {stroke} bind def\n"
             "/sh {show} bind def\n"
             "/slc {setlinecap} bind def\n"
             "/slj {setlinejoin} bind def\n"
             "/slw {setlinewidth} bind def\n"
             "/srgb {setrgbcolor} bind def\n"
             "/rot {rotate} bind def\n"
             "/sc {scale} bind def\n"
             "/sd {setdash} bind def\n"
             "/ff {findfont} bind def\n"
             "/sf {setfont} bind def\n"
             "/scf {scalefont} bind def\n"
             "/sw {stringwidth pop} bind def\n"
             "/tr {translate} bind def\n"

             "/JR {\n"
             " neg 0\n"
             " rmoveto\n"
             "} bind def\n"

             "/JC {\n"
             " 2 div neg 0\n"
             " rmoveto\n"
             "} bind def\n"
  
             "\n/ellipsedict 8 dict def\n"
             "ellipsedict /mtrx matrix put\n"
             "/ellipse\n"
             "{ ellipsedict begin\n"
             "   /endangle exch def\n"
             "   /startangle exch def\n"
             "   /yrad exch def\n"
             "   /xrad exch def\n"
             "   /y exch def\n"
             "   /x exch def"
             "   /savematrix mtrx currentmatrix def\n"
             "   x y tr xrad yrad sc\n"
             "   0 0 1 startangle endangle arc\n"
             "   savematrix setmatrix\n"
             "   end\n"
             "} def\n\n"
    ); 
  
    fprintf(psout,
          "[ /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef\n"
          "/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef\n"
          "/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef\n"
          "/.notdef /.notdef /space /exclam /quotedbl /numbersign /dollar /percent /ampersand /quoteright\n"
          "/parenleft /parenright /asterisk /plus /comma /hyphen /period /slash /zero /one\n"
          "/two /three /four /five /six /seven /eight /nine /colon /semicolon\n"          "/less /equal /greater /question /at /A /B /C /D /E\n"
          "/F /G /H /I /J /K /L /M /N /O\n"
          "/P /Q /R /S /T /U /V /W /X /Y\n"
          "/Z /bracketleft /backslash /bracketright /asciicircum /underscore /quoteleft /a /b /c\n"
          "/d /e /f /g /h /i /j /k /l /m\n"
          "/n /o /p /q /r /s /t /u /v /w\n"
          "/x /y /z /braceleft /bar /braceright /asciitilde /.notdef /.notdef /.notdef\n"
          "/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef\n"
          "/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef\n"
          "/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef\n"
          "/space /exclamdown /cent /sterling /currency /yen /brokenbar /section /dieresis /copyright\n"
          "/ordfeminine /guillemotleft /logicalnot /hyphen /registered /macron /degree /plusminus /twosuperior /threesuperior\n"
          "/acute /mu /paragraph /periodcentered /cedilla /onesuperior /ordmasculine /guillemotright /onequarter /onehalf\n"
          "/threequarters /questiondown /Agrave /Aacute /Acircumflex /Atilde /Adieresis /Aring /AE /Ccedilla\n"
          "/Egrave /Eacute /Ecircumflex /Edieresis /Igrave /Iacute /Icircumflex /Idieresis /Eth /Ntilde\n"
          "/Ograve /Oacute /Ocircumflex /Otilde /Odieresis /multiply /Oslash /Ugrave /Uacute /Ucircumflex\n"
          "/Udieresis /Yacute /Thorn /germandbls /agrave /aacute /acircumflex /atilde /adieresis /aring\n"
          "/ae /ccedilla /egrave /eacute /ecircumflex /edieresis /igrave /iacute /icircumflex /idieresis\n"
          "/eth /ntilde /ograve /oacute /ocircumflex /otilde /odieresis /divide /oslash /ugrave\n"
          "/uacute /ucircumflex /udieresis /yacute /thorn /ydieresis] /isolatin1encoding exch def\n");
 
    ps_reencode_font(psout, "Times-Roman");
    ps_reencode_font(psout, "Times-Italic");
    ps_reencode_font(psout, "Times-Bold");
    ps_reencode_font(psout, "Times-BoldItalic");
    ps_reencode_font(psout, "AvantGarde-Book");
    ps_reencode_font(psout, "AvantGarde-BookOblique");
    ps_reencode_font(psout, "AvantGarde-Demi");
    ps_reencode_font(psout, "AvantGarde-DemiOblique");
    ps_reencode_font(psout, "Bookman-Light");
    ps_reencode_font(psout, "Bookman-LightItalic");
    ps_reencode_font(psout, "Bookman-Demi");
    ps_reencode_font(psout, "Bookman-DemiItalic");
    ps_reencode_font(psout, "Courier");
    ps_reencode_font(psout, "Courier-Oblique");
    ps_reencode_font(psout, "Courier-Bold");
    ps_reencode_font(psout, "Courier-BoldOblique");
    ps_reencode_font(psout, "Helvetica");
    ps_reencode_font(psout, "Helvetica-Oblique");
    ps_reencode_font(psout, "Helvetica-Bold");
    ps_reencode_font(psout, "Helvetica-BoldOblique");
    ps_reencode_font(psout, "Helvetica-Narrow");
    ps_reencode_font(psout, "Helvetica-Narrow-Oblique");
    ps_reencode_font(psout, "Helvetica-Narrow-Bold");
    ps_reencode_font(psout, "Helvetica-Narrow-BoldOblique");
    ps_reencode_font(psout, "NewCenturySchoolbook-Roman");
    ps_reencode_font(psout, "NewCenturySchoolbook-Italic");
    ps_reencode_font(psout, "NewCenturySchoolbook-Bold");
    ps_reencode_font(psout, "NewCenturySchoolbook-BoldItalic");
    ps_reencode_font(psout, "Palatino-Roman");
    ps_reencode_font(psout, "Palatino-Italic");
    ps_reencode_font(psout, "Palatino-Bold");
    ps_reencode_font(psout, "Palatino-BoldItalic");
    ps_reencode_font(psout, "Symbol");
    ps_reencode_font(psout, "ZapfChancery-MediumItalic");
    ps_reencode_font(psout, "ZapfDingbats");
   
    if(ps->orientation == GTK_PLOT_PORTRAIT)
             fprintf(psout, "%d %d translate\n"
                            "%g %g scale\n",
                            0, ps->page_height,
                            ps->scalex, -ps->scaley);

    if(ps->orientation == GTK_PLOT_LANDSCAPE)
             fprintf(psout, "%g %g scale\n"
                            "-90 rotate \n",
                            ps->scalex, -ps->scaley);

    fprintf(psout,"%%%%EndProlog\n\n\n");

    return TRUE;
}

static void ps_reencode_font(FILE *file, char *fontname)
{
  if (!strcmp(fontname, "Symbol"))
    fprintf(file,
            "/%s-latin1\n"
            "    /%s findfont\n"
            "definefont pop\n", fontname, fontname);
  else
    fprintf(file,
            "/%s-latin1\n"
            "    /%s findfont\n"
            "    dup length dict begin\n"
            "   {1 index /FID ne {def} {pop pop} ifelse} forall\n"
            "   /Encoding isolatin1encoding def\n"
            "    currentdict end\n"
            "definefont pop\n", fontname, fontname);
}

/* ............................................................. */

static int ps_print_plots (BOXPLOT **plotgrp, int nplots, 
			   int winwidth, int winheight,
			   double gmax, double gmin,
			   GtkStyle *style)
{
    GtkPlotPS *ps;
    int i;

    ps = GTK_PLOT_PS(gtk_plot_ps_new ("foo.eps", 
				      GTK_PLOT_PORTRAIT, 
				      TRUE, /* epsflag */
				      GTK_PLOT_LETTER, 
				      1.0, 1.0)
		     );

    ps->pc.width = winwidth;
    ps->pc.height = winheight;

    psinit (ps);
    gtk_psfont_init();

    for (i=0; i<nplots; i++)
	gtk_area_boxplot (plotgrp[i], gmax, gmin, NULL, 
			  style, &ps->pc, winheight);
    
    gtk_boxplot_yscale(NULL, style, &ps->pc, gmax, gmin, winheight);

    psleave(ps);
    gtk_object_destroy(GTK_OBJECT(ps));
    gtk_psfont_unref();

    return 0;
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
			  boxwin->style, NULL, winheight);
    
    gtk_boxplot_yscale(area, boxwin->style, NULL, gmax, gmin, winheight);

    /* just testing -- needs to be put on a menu or something */
    gdk_color_white (gdk_colormap_get_system(), &white);
    gdk_color_black (gdk_colormap_get_system(), &black);

    ps_print_plots (plotgrp, nplots, winwidth, winheight, 
		    gmax, gmin, boxwin->style);

    free_plotgrp(plotgrp);

    return 0;
}




