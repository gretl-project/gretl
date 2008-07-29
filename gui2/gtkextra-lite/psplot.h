/* postscript driver, based on gtkplot code,
 * Copyright 1999-2001  Adrian E. Feiguin <feiguin@ifir.edu.ar>
 * Slimmed down for use with gretl by Allin Cottrell <cottrell@wfu.edu>
 *
 * Some few lines of code borrowed from
 * DiaCanvas -- a technical canvas widget
 * Copyright (C) 1999 Arjan Molenaar
 * Dia -- an diagram creation/manipulation program
 * Copyright (C) 1998 Alexander Larsson
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef PS_PLOT_H_
#define PS_PLOT_H_

enum{
    GTK_PLOT_LETTER,
    GTK_PLOT_A4
};

enum{
    GTK_PLOT_PORTRAIT,
    GTK_PLOT_LANDSCAPE	
};

enum {
    GTK_PLOT_BORDER_NONE,
    GTK_PLOT_BORDER_LINE,
    GTK_PLOT_BORDER_SHADOW
};

typedef struct _PSPlot PSPlot;
typedef struct _PlotPoint PlotPoint;

struct _PlotPoint {
    gdouble x, y;
};

void ps_plot_init (PSPlot *ps);

void ps_plot_finalize (PSPlot *ps); 

PSPlot *ps_plot_new (const gchar *fname,
		     FILE *fp,
		     gint orientation,
		     gint epsflag,
		     gint page_size,
		     gdouble scalex,
		     gdouble scaley);

void ps_plot_set_page_size (PSPlot *ps, gdouble width, gdouble height);

void ps_plot_draw_line (PSPlot *ps, gdouble x0, gdouble y0, 
			gdouble xf, gdouble yf);

void ps_plot_draw_circle (PSPlot *ps, gboolean filled, 
			  gdouble x, gdouble y, gdouble size);

void ps_plot_draw_polygon (PSPlot *ps, gboolean filled, PlotPoint *points, 
			   gint numpoints);

void ps_plot_draw_rectangle (PSPlot *ps, gboolean filled, gdouble x, gdouble y, 
			     gdouble width, gdouble height);

void ps_plot_draw_string (PSPlot *ps,
			  gint x, gint y,
			  gint angle,
			  const GdkColor *fg,
			  const GdkColor *bg,
			  gboolean transparent,
			  gint border,
			  gint border_space,
			  gint border_width,
			  gint shadow_width,
			  const gchar *font,
			  gint font_height,
			  GtkJustification justification,
			  const gchar *text);

#endif
