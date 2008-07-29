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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <gtk/gtk.h>

#include "gtkpsfontpango.h"
#include "psplot.h"

#define PLOT_LETTER_W 	612   /* Width and Height in ps points */
#define PLOT_LETTER_H 	792

#define PLOT_A4_W       595
#define PLOT_A4_H       842

struct _PSPlot
{
   FILE *psfile;
   gchar *psname;
 
   gint orientation;
   gint epsflag;

   gint page_size;
   gint page_width;
   gint page_height;

   gdouble scalex, scaley;

   gboolean gsaved;
};

PSPlot *ps_plot_new (const gchar *fname,
		     FILE *fp,
		     gint orientation,
		     gint epsflag,
		     gint page_size,
		     gdouble epswidth,
		     gdouble epsheight,
		     gdouble scalex,
		     gdouble scaley)
{
    PSPlot *ps = malloc(sizeof *ps);

    if (ps == NULL) {
	return NULL;
    }

    ps->psname = g_strdup(fname);
    ps->psfile = fp;
    ps->gsaved = FALSE;
    ps->orientation = orientation;
    ps->epsflag = epsflag;
    ps->page_size = page_size;
    ps->scalex = scalex;
    ps->scaley = scaley;

    if (epsflag) {
	ps->page_width = epswidth;
	ps->page_height = epsheight;
    } else if (page_size == PS_PLOT_A4) {
        ps->page_width = PLOT_A4_W;
        ps->page_height = PLOT_A4_H;
    } else {	
        ps->page_width = PLOT_LETTER_W;
        ps->page_height = PLOT_LETTER_H;
    }

    return ps;
}

static void pssetcolor (PSPlot *ps, const GdkColor *color)
{
    fprintf(ps->psfile, "%g %g %g setrgbcolor\n",
	    (gdouble) color->red / 65535.0,
	    (gdouble) color->green / 65535.0,
	    (gdouble) color->blue / 65535.0);
}

static void pssetlineattr (PSPlot *ps, 
			   gfloat line_width,
			   GdkLineStyle line_style,
			   GdkCapStyle cap_style,
			   GdkJoinStyle join_style)
{
    fprintf(ps->psfile, "%g slw\n", line_width);
    fprintf(ps->psfile, "%d slc\n", abs(cap_style - 1));
    fprintf(ps->psfile, "%d slj\n", join_style);

    if (line_style == 0) {
	fputs("[] 0 sd\n", ps->psfile);  /* solid line */
    }
}

static void pssetdash (PSPlot *ps,
		       gdouble offset, 
		       gdouble *values,
		       gint num_values)
{
    FILE *psout = ps->psfile;

    switch (num_values) {
    case 0:
        fputs("[] 0 sd\n", psout);
        break;
    case 2:
        fprintf(psout, "[%g %g] %g sd\n", values[0], values[1], offset);
        break;
    case 4:
        fprintf(psout, "[%g %g %g %g] %g sd\n", values[0], values[1],
		values[2], values[3], 
		offset);
        break;
    case 6:
        fprintf(psout, "[%g %g %g %g %g %g] %g sd\n",  values[0], values[1],
		values[2], values[3], 
		values[4], values[5], 
		offset);
        break;
    default:
        break;
    }
}

void ps_plot_draw_line (PSPlot *ps, gdouble x0, gdouble y0, 
			gdouble xf, gdouble yf)
{
    fprintf(ps->psfile, "%g %g m\n", x0, y0);
    fprintf(ps->psfile, "%g %g l\n", xf, yf);
    fputs("s\n", ps->psfile);
}

void ps_plot_draw_circle (PSPlot *ps, gboolean filled, gdouble x, 
			  gdouble y, gdouble size)
{
    fprintf(ps->psfile,"n %g %g %g %g 0 360 ellipse\n", 
	    x, y, size / 2., size / 2.);

    if (filled) {
	fputs("f\n", ps->psfile);
    }

    fputs("s\n", ps->psfile);
}

void ps_plot_draw_polygon (PSPlot *ps, gboolean filled, PlotPoint *points, 
			   gint numpoints)
{
    FILE *psout = ps->psfile;
    gint i;

    fputs("n\n", psout);
    fprintf(psout, "%g %g m\n", points[0].x, points[0].y);

    for (i=1; i<numpoints; i++) {
	fprintf(psout,"%g %g l\n", points[i].x, points[i].y);
    }

    if (filled) {
	fputs("f\n", psout);
    } else {
	fputs("cp\n", psout);
    }

    fputs("s\n", psout);
}

void ps_plot_draw_rectangle (PSPlot *ps, gboolean filled, gdouble x, gdouble y, 
			     gdouble width, gdouble height)
{
    PlotPoint point[4];

    point[0].x = x;
    point[0].y = y;
    point[1].x = x + width;
    point[1].y = y;
    point[2].x = x + width;
    point[2].y = y + height;
    point[3].x = x;
    point[3].y = y + height;

    ps_plot_draw_polygon(ps, filled, point, 4);
}

static void pssetfont (PSPlot *ps, GtkPSFont *psfont, gint height)
{
    FILE *psout = ps->psfile;

    if (psfont->i18n_latinfamily && psfont->vertical) {
	fprintf(psout,
		"/%s ff [0 1 -1 0 0 0.3] makefont [%d 0 0 %d 0 0] makefont sf\n",
		psfont->psname, height, height);
    } else {
	fprintf(psout, "/%s-latin1 ff %g scf sf\n", psfont->psname, (double) height);
    }
}

static void psgsave (PSPlot *ps)
{
    fputs("gsave\n", ps->psfile);
    ps->gsaved = TRUE;
}

static void psgrestore (PSPlot *ps)
{
    if (ps->gsaved) {
	fputs("grestore\n", ps->psfile);
	ps->gsaved = FALSE;
    }
}

static void
psoutputstring (PSPlot *ps,
		GtkPSFont *psfont,
		GtkPSFont *latin_psfont,
		gint height,
		const GdkWChar *wstring,
		const gchar *addstring)
{
    const GdkWChar *p;
    GdkWChar wcs[2];
    gint curcode = 0, code; /* 0..neither 1..latin 2..i18n */
    gchar begin[] = { 0, '(', '<' };
    gchar end[] = { 0, ')', '>' };
    GtkPSFont *fonts[] = { NULL, latin_psfont, psfont };
    FILE *out = ps->psfile;
    gchar *mbs, *c;

    p = wstring;
  
    if (psfont->i18n_latinfamily) {
	while (*p) {
	    code = (*p <= 0x7f) ? 1 : 2;
	    if (curcode && curcode != code)
		fprintf(out, "%c %s\n", end[curcode], addstring);
	    if (curcode != code) {
		pssetfont(ps, fonts[code], height);
		fputc(begin[code], out);
	    }
      
	    curcode = code;

	    wcs[0] = *p++;
	    wcs[1] = 0;
	    c = mbs = gdk_wcstombs(wcs);
      
	    if (code == 2) {
		while (*c) {
		    fprintf(out, "%02x", (unsigned char)(*c++));
		}
	    } else {
		if (*c == '(' || *c == ')') {
		    fputc('\\', out);
		}
		fputc(*c, out);
	    }
      
	    g_free(mbs);
	}
    } else {
	c = mbs = gdk_wcstombs(wstring);
    
	while (*c) {
	    if (curcode == 0) {
		pssetfont(ps, psfont, height);
		fputc(begin[1], out);
		curcode = 1;
	    }

	    if (*c == '(' || *c == ')') {
		fputc('\\', out);
	    }
	    fputc(*c++, out);
	}

	g_free(mbs);
    }

    if (curcode) {
	fprintf(out, "%c %s\n", end[curcode], addstring);
    }
}

static void maybe_nullify_latin (GdkFont **latin_font)
{
    if (*latin_font) {
	gdk_font_unref(*latin_font);
	*latin_font = NULL;
    }
}

static void
plot_text_get_size(const gchar *text, gint angle, 
		   const gchar* text_font, 
		   gint text_height,   
		   gint *width, gint *height,
		   gint *ascent, gint *descent)
{
    GdkFont *font = NULL, *latin_font = NULL;
    GtkPSFont *psfont, *base_psfont, *latin_psfont;
    gint old_width, old_height;
    gboolean italic = FALSE, bold = FALSE;
    gint fontsize;
    gint x, y, y0, w;
    GList *family;
    gint numf;
    GdkWChar *aux, *wtext, *lastchar = NULL;
    gint i, a, d;

    gtk_psfont_get_families(&family, &numf);
    base_psfont = psfont = gtk_psfont_get_by_name(text_font);
    font = gtk_psfont_get_gdkfont(psfont, text_height);
    old_width = gdk_string_width (font, text);
    old_height = font->ascent + font->descent;

    if (psfont->i18n_latinfamily) {
	latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily, italic,
						bold);
	latin_font = gtk_psfont_get_gdkfont(latin_psfont, text_height);
    } else {
	latin_font = NULL;
	latin_psfont = NULL;
    }

    italic = psfont->italic;
    bold = psfont->bold;
    fontsize = text_height;
  
    x = 0;
    y0 = y = font->ascent;
    old_width = 0;

    *ascent = font->ascent;
    *descent = font->descent;

    i = strlen(text) + 2;
    aux = wtext = g_malloc0(i * sizeof(GdkWChar));
    gdk_mbstowcs(wtext, text, i - 1);

    while (aux && *aux != '\0' && *aux != '\n') {
	if (*aux == '\\'){
	    aux++;
	    switch (*aux) {
	    case '0': case '1': case '2': case '3':
	    case '4': case '5': case '6': case '7': case '9':
		psfont = gtk_psfont_get_by_family((gchar *) g_list_nth_data(family, *aux-'0'), 
						  italic, bold);
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		aux++;
		break;
	    case '8': case 'g':
		psfont = gtk_psfont_get_by_family("Symbol", italic, bold);
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		aux++;
		break;
	    case 'B':
		bold = TRUE;
		psfont = gtk_psfont_get_by_family(psfont->family, italic, bold);
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		maybe_nullify_latin(&latin_font);
		if (psfont->i18n_latinfamily) {
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily,
							    italic, bold);
		    latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
		}
		aux++;
		break;
	    case 'i':
		italic = TRUE;
		psfont = gtk_psfont_get_by_family(psfont->family, italic, bold);
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		maybe_nullify_latin(&latin_font);
		if (psfont->i18n_latinfamily) {
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily,
							    italic, bold);
		    latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
		}
		aux++;
		break;
	    case 'S': case '^':
		fontsize = (int) ((gdouble)fontsize * 0.6 + 0.5);
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		maybe_nullify_latin(&latin_font);
		if (psfont->i18n_latinfamily) {
		    latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
		}
		y -= font->ascent;
		aux++;
		break;
	    case 's': case '_':
		fontsize = (int) ((gdouble)fontsize * 0.6 + 0.5);
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		maybe_nullify_latin(&latin_font);
		if (psfont->i18n_latinfamily) {
		    latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
		}
		y += font->descent;
		aux++;
		break;
	    case '+':
		fontsize += 3;
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		maybe_nullify_latin(&latin_font);
		if (psfont->i18n_latinfamily) {
		    latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
		}
		aux++;
		break;
	    case '-':
		fontsize -= 3;
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, fontsize);
		maybe_nullify_latin(&latin_font);
		if (psfont->i18n_latinfamily) {
		    latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
		}
		aux++;
		break;
	    case 'N':
		psfont = base_psfont;
		gdk_font_unref(font);
		font = gtk_psfont_get_gdkfont(psfont, text_height);
		y = y0;
		italic = psfont->italic;
		bold = psfont->bold;
		fontsize = text_height;
		maybe_nullify_latin(&latin_font);
		if (psfont->i18n_latinfamily) {
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily,
							    italic, bold);
		    latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
		}
		aux++;
		break;
	    case 'b':
		if (lastchar) {
		    gtk_psfont_get_char_size(psfont, font, latin_font, *lastchar, &w,
					     NULL, NULL);
		    x -= w;
	     
		    if (lastchar == wtext) {
			lastchar = NULL;
		    } else {
			lastchar--;
		    }
		} else {
		    gtk_psfont_get_char_size(psfont, font, latin_font, 'X', &w, NULL, NULL);
		    x -= w;
		}
		aux++;
		break;
	    default:
		if (aux && *aux != '\0' && *aux !='\n') {
		    gtk_psfont_get_char_size(psfont, font, latin_font, *aux, &w, &a, &d);
		    x += w;
		    lastchar = aux;
		    aux++;
		}
		break;
	    }
	} else {
	    if (aux && *aux != '\0' && *aux != '\n') {
		gtk_psfont_get_char_size(psfont, font, latin_font, *aux, &w, &a, &d);
		x += w;
		lastchar = aux;
		aux++;
		if (x > old_width) old_width = x;
		if (y + d - y0 > *descent) *descent = y + d - y0;
		if (y0 - y + a > *ascent) *ascent = y0 - y + a;
	    }
	}
    }

    old_height = *ascent + *descent;
    *width = old_width;
    *height = old_height;
    if (angle == 90 || angle == 270) {
	*width = old_height;
	*height = old_width;
    }

    g_free(wtext);
    gdk_font_unref(font);
}

static void
ps_plot_text_get_area(const gchar *text, gint angle, GtkJustification just, 
		      const gchar *font, gint font_height,
		      gint *x, gint *y,
		      gint *width, gint *height)
{
    gint ascent, descent;

    if (text == NULL) return;

    plot_text_get_size(text, angle, font, font_height, 
		       width, height, &ascent, &descent);

    *x = 0;
    *y = 0;

    switch (just) {
    case GTK_JUSTIFY_LEFT:
	switch (angle){
        case 0:
            *y -= ascent;
            break;
        case 90:
            *y -= *height;
            *x -= ascent;
            break;
        case 180:
            *x -= *width;
            *y -= descent;
            break;
        case 270:
            *x -= descent;
            break;
	}
	break;
    case GTK_JUSTIFY_RIGHT:
	switch (angle) {
        case 0:
            *x -= *width;
            *y -= ascent;
            break;
        case 90:
            *x -= ascent;
            break;
        case 180:
            *y -= descent; 
            break;
        case 270:
            *y -= *height;
            *x -= descent; 
            break;
	}
	break;
    case GTK_JUSTIFY_CENTER:
    default:
	switch (angle) {
        case 0:
            *x -= *width / 2.;
            *y -= ascent;
	    break;
        case 90:
            *x -= ascent;
            *y -= *height / 2.;
	    break;
        case 180:
            *x -= *width / 2.;
            *y -= descent;
            break;
        case 270:
            *x -= descent;
            *y -= *height / 2.;
            break;
	}
    }
}

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
			  const gchar *text)
{
    gchar *currfont;
    const gchar *c;
    GtkPSFont *psfont, *base_psfont, *latin_psfont = NULL;
    gint curcnt = 0, offset = 0;
    gint numf;
    gdouble scale;
    gboolean italic, bold;
    gboolean special = FALSE;
    GList *family;
    FILE *psout;
    gint twidth, theight, tdescent, tascent;
    gint tx, ty, width, height; 
    gint i;
    GdkWChar *curstr, *wtext, *lastchar = NULL, bkspchar[2], *aux, *xaux;
    gchar num[4];

    if (text == NULL || *text == '\0') {
	return;
    }

    psout = ps->psfile;

    gtk_psfont_get_families(&family, &numf);
    base_psfont = psfont = gtk_psfont_get_by_name(font);
    italic = psfont->italic;
    bold = psfont->bold;

    currfont = psfont->psname;

    if (psfont->i18n_latinfamily) {
	latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily, italic,
						bold);
    }

    ps_plot_text_get_area(text, angle, justification, font, font_height,
			  &tx, &ty, &width, &height);

    tx += x;
    ty += y;

    if (!transparent) {
	pssetcolor(ps, bg);
	ps_plot_draw_rectangle(ps,
			       TRUE,
			       tx - border_space, ty - border_space,
			       width + 2*border_space, height + 2*border_space);
    }
    /* border */

    pssetcolor(ps, fg);
    pssetdash(ps, 0, NULL, 0);
    pssetlineattr(ps, border_width, 0, 0, 0);
 
    switch (border){
    case PS_PLOT_BORDER_SHADOW:
	ps_plot_draw_rectangle(ps, TRUE, 
			       tx - border_space + shadow_width,
			       ty + height + border_space, 
			       width + 2 * border_space, shadow_width);
	ps_plot_draw_rectangle(ps, TRUE, 
			       tx + width + border_space, 
			       ty - border_space + shadow_width, 
			       shadow_width, height + 2 * border_space);
    case PS_PLOT_BORDER_LINE: 
	ps_plot_draw_rectangle(ps, FALSE, 
			       tx - border_space, ty - border_space, 
			       width + 2*border_space, height + 2*border_space);
    case PS_PLOT_BORDER_NONE:
    default:
        break;
    }

    plot_text_get_size(text, angle, psfont->psname, font_height, 
		       &twidth, &theight, &tascent, &tdescent);

    if (angle == 90 || angle == 270) {
	angle = 360 - angle;
    }

    psgsave(ps);
    fprintf(psout, "%d %d translate\n", x, y);
    fprintf(psout, "%d rotate\n", angle);

    fprintf(psout, "0 0 m\n");
    fprintf(psout, "1 -1 sc\n");

    if (psfont->i18n_latinfamily) {
	special = TRUE;
    }

    c = text;
    while (c && *c != '\0' && *c != '\n') {
	if(*c == '\\'){
	    c++;
	    switch (*c) {
	    case '0': case '1': case '2': case '3':
	    case '4': case '5': case '6': case '7': case '9':
	    case '8': case'g': case 'B': case 'b': case 'x': case 'N':
	    case 's': case 'S': case 'i': case '-': case '+': case '^':
		special = TRUE;
		break;
	    default:
		break;
	    }
	} else {
	    c++;
	}
    }

    if (special) {
	switch (justification) {
	case GTK_JUSTIFY_LEFT:
	    break;
	case GTK_JUSTIFY_RIGHT:
	    if (angle == 0 || angle == 180) {
		fprintf(psout, "%d JR\n", twidth);
	    } else {
		fprintf(psout, "%d JR\n", theight);
	    }
	    break;
	case GTK_JUSTIFY_CENTER:
	default:
	    if (angle == 0 || angle == 180) {
		fprintf(psout, "%d JC\n", twidth);
	    } else {
		fprintf(psout, "%d JC\n", theight);
	    }
	    break;
	}
    } else {
	pssetfont(ps, psfont, font_height);
    
	switch (justification) {
	case GTK_JUSTIFY_LEFT:
	    break;
	case GTK_JUSTIFY_RIGHT:
	    fprintf(psout, "(%s) sw JR\n", text);
	    break;
	case GTK_JUSTIFY_CENTER:
	default:
	    fprintf(psout, "(%s) sw JC\n", text);
	    break;
	}
	fprintf(psout, "(%s) show\n", text);

	psgrestore(ps);  
	fputs("n\n", psout);  
	return;
    }

    i = strlen(text) + 2;
    curstr = g_malloc0(i * sizeof(GdkWChar));
    aux = wtext = g_malloc0(i * sizeof(GdkWChar));
    gdk_mbstowcs(wtext, text, i - 1);

    scale = font_height;
    curcnt = 0;

    while (aux && *aux != '\0' && *aux != '\n') {
	if (*aux == '\\'){
	    aux++;
	    switch(*aux){
	    case '0': case '1': case '2': case '3':
	    case '4': case '5': case '6': case '7': case '9':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		psfont = gtk_psfont_get_by_family((gchar *)g_list_nth_data(family, *aux-'0'), 
						  italic, bold);
		aux++;
		break;
	    case '8': case 'g':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		psfont = gtk_psfont_get_by_family("Symbol", italic, bold);
		aux++;
		break;
	    case 'B':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		bold = TRUE;
		psfont = gtk_psfont_get_by_family(psfont->family, italic, bold);
		if (psfont->i18n_latinfamily) {
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily, italic, bold);
		}
		aux++;
		break;
	    case 'x':
		xaux = aux + 1;
		for (i=0; i<3; i++) {
		    if (xaux[i] >= '0' && xaux[i] <= '9') {
			num[i] = xaux[i];
		    } else {
			break;
		    }
		}
		if (i < 3){
		    aux++;
		    break;
		}
		num[3] = '\0';
		  
		i = atoi(num);
		g_snprintf(num, 4, "%o", i % (64 * 8));

		curstr[curcnt++] = '\\';
		i = 0;
		while (num[i]) {
		    curstr[curcnt++] = num[i++];
		}
		  
		aux += 4;
		break;
	    case 'i':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		italic = TRUE;
		psfont = gtk_psfont_get_by_family(psfont->family, italic, bold);
		if (psfont->i18n_latinfamily)
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily, italic, bold);
		aux++;
		break;
	    case 's': case '_':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		scale = 0.6 * font_height;
		offset -= (gint)scale / 2;
		fprintf(psout, "0 %d rmoveto\n", -((gint)scale / 2));
		aux++;
		break;
	    case 'S': case '^':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		scale = 0.6 * font_height;
		offset += 0.5 * font_height;
		fprintf(psout, "0 %d rmoveto\n", (gint)(0.5*font_height));
		aux++;
		break;
	    case 'N':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		psfont = base_psfont;
		italic = psfont->italic;
		bold = psfont->bold;
		if (psfont->i18n_latinfamily) {
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily,
							    italic, bold);
		}
		scale = font_height;
		fprintf(psout, "0 %d rmoveto\n", -offset);
		offset = 0;
		aux++;
		break;
	    case 'b':
		curstr[curcnt] = '\0';
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		bkspchar[1] = 0;
		if (lastchar) {
		    bkspchar[0] = *lastchar;
		    lastchar--;
		} else {
		    bkspchar[0] = 'X';
		    lastchar = NULL;
		}
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       bkspchar,
			       "stringwidth pop 0 exch neg exch rmoveto");
		aux++;
		break;
	    case '-':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		scale -= 3;
		if (scale < 6) {
		    scale = 6;
		}
		aux++;
		break;
	    case '+':
		curstr[curcnt] = 0;
		psoutputstring(ps, psfont, latin_psfont, (gint)scale,
			       curstr, "show");
		curcnt = 0;
		scale += 3;
		aux++;
		break;
	    default:
		if(aux && *aux != '\0' && *aux != '\n'){
                    curstr[curcnt++] = *aux;
                    aux++;
		}
		break;
	    }
	} else if (aux && *aux != '\0' && *aux != '\n') {
	    curstr[curcnt++] = *aux;
	    lastchar = aux;
	    aux++;
	}
    }

    curstr[curcnt] = 0;
    psoutputstring(ps, psfont, latin_psfont, (gint) scale, curstr, "show");

    psgrestore(ps);  
    fputs("n\n", psout);  

    g_free(wtext);
    g_free(curstr);
}

static void ps_reencode_font (FILE *file, char *fontname)
{
    /* Don't reencode the Symbol font, as it doesn't work in latin1 encoding.
     * Instead, just define Symbol-latin1 to be the same as Symbol. */
    if (!strcmp(fontname, "Symbol")) {
	fprintf(file,
		"/%s-latin1\n"
		"    /%s findfont\n"
		"definefont pop\n", fontname, fontname);
    } else {
	fprintf(file,
		"/%s-latin1\n"
		"    /%s findfont\n"
		"    dup length dict begin\n"
		"   {1 index /FID ne {def} {pop pop} ifelse} forall\n"
		"   /Encoding isolatin1encoding def\n"
		"    currentdict end\n"
		"definefont pop\n", fontname, fontname);
    }
}

void ps_plot_init (PSPlot *ps)
{
    time_t now = time(NULL);
    FILE *psout = ps->psfile;

    if (ps->epsflag) {
	fputs("%%!PS-Adobe-2.0 PCF-2.0\n", psout);
    } else {
	fputs("%%!PS-Adobe-2.0\n", psout);
    }

    fprintf(psout,
	    "%%%%Title: %s\n"
	    "%%%%Creator: %s v%s Copyright (c) 1999 Adrian E. Feiguin\n"
	    "%%%%CreationDate: %s"
	    "%%%%Magnification: 1.0000\n",
	    ps->psname,
	    "GtkPlot", "3.x",
	    ctime(&now));

    if (ps->orientation == PS_PLOT_PORTRAIT) {
	fputs("%%%%Orientation: Portrait\n", psout);
    } else {
	fputs("%%%%Orientation: Landscape\n", psout);
    }

    if (ps->epsflag) {
	fprintf(psout,
		"%%%%BoundingBox: 0 0 %d %d\n"
		"%%%%Pages: 1\n"
		"%%%%EndComments\n",
		ps->page_width,
		ps->page_height);
    }

    fputs ("/cp {closepath} bind def\n"
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
	   "} def\n\n", psout);
  
    fputs("[ /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef\n"
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
	  "/uacute /ucircumflex /udieresis /yacute /thorn /ydieresis] /isolatin1encoding exch def\n", psout);
 
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
   
    if (ps->orientation == PS_PLOT_PORTRAIT) {
	fprintf(psout, "%d %d translate\n"
		"%g %g scale\n",
		0, ps->page_height,
		ps->scalex, -ps->scaley);
    } else if (ps->orientation == PS_PLOT_LANDSCAPE) {
	fprintf(psout, "%g %g scale\n"
		"-90 rotate \n",
		ps->scalex, -ps->scaley);
    }

    fputs("%%%%EndProlog\n\n\n", psout);

    gtk_psfont_init();
}

void ps_plot_finalize (PSPlot *ps)
{
    FILE *fp = ps->psfile;

    fputs("showpage\n", fp);
    fputs("%%%%Trailer\n", fp);
    fputs("%%%%EOF\n", fp);

    gtk_psfont_unref();

    g_free(ps->psname);
    free(ps);
}
