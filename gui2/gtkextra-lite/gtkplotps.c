/* gtkplotps - postscript driver
 * Copyright 1999-2001  Adrian E. Feiguin <feiguin@ifir.edu.ar>
 *
 * Some few lines of code borrowed from
 * DiaCanvas -- a technical canvas widget
 * Copyright (C) 1999 Arjan Molenaar
 * Dia -- an diagram creation/manipulation program
 * Copyright (C) 1998 Alexander Larsson
 * ISO Latin encoding by
 * Przemek Klosowski
 * przemek@rrdbartok.nist.gov
 * (borrowed from XMGR)
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

#include "gretl_intl.h"

#include "gtkplot-lite.h"
#include "gtkpsfontpango.h"
#include "gtkplotpc.h"
#include "gtkplotps.h"

/* I've made these public for use with gretl -- Allin Cottrell */
void psleave				(GtkPlotPC *pc);
gboolean psinit				(GtkPlotPC *pc); 

static void
gtk_plot_text_get_size(const gchar *text, gint angle, 
                       const gchar* text_font, 
                       gint text_height,   
                       gint *width, gint *height,
                       gint *ascent, gint *descent);
static void
gtk_plot_text_get_area(const gchar *text, gint angle, GtkJustification just, 
                       const gchar *font, gint font_height,
                       gint *x, gint *y,
                       gint *width, gint *height);


static void gtk_plot_ps_class_init 		(GtkPlotPSClass *klass);
static void gtk_plot_ps_init 			(GtkPlotPS *ps);
static void gtk_plot_ps_destroy 		(GtkObject *object);

/*********************************************************************/
/* Postscript specific functions */

static void pssetviewport			(GtkPlotPC *pc, 
						 gdouble w, gdouble h); 
static void psgsave				(GtkPlotPC *pc);
static void psgrestore				(GtkPlotPC *pc);
static void psdrawlines				(GtkPlotPC *pc,
						 GtkPlotPoint *points, 
						 gint numpoints);
static void psdrawpoint				(GtkPlotPC *pc, 
                				 gdouble x, gdouble y); 
static void psdrawline				(GtkPlotPC *pc,
						 gdouble x0, gdouble y0, 
						 gdouble xf, gdouble yf);
static void psdrawpolygon			(GtkPlotPC *pc,
						 gboolean filled,
						 GtkPlotPoint *points, 
						 gint numpoints); 
static void psdrawrectangle			(GtkPlotPC *pc, 
						 gboolean filled, 
                				 gdouble x, gdouble y, 
						 gdouble width, gdouble height);
static void psdrawcircle			(GtkPlotPC *pc,
						 gboolean filled,
                                                 gdouble x, gdouble y, 
						 gdouble size);
static void psdrawellipse			(GtkPlotPC *pc, 
              					 gboolean filled,
						 gdouble x, gdouble y, 
						 gdouble width, gdouble height); 
static void pssetcolor				(GtkPlotPC *pc, 
						 const GdkColor *color); 
static void pssetlineattr			(GtkPlotPC *pc, 
                                                 gfloat line_width,
                                                 GdkLineStyle line_style,
                                                 GdkCapStyle cap_style,
                                                 GdkJoinStyle join_style);
static void psdrawstring			(GtkPlotPC *pc,
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
                                                 gint height,
                                                 GtkJustification just,
                                                 const gchar *text);
static void pssetfont				(GtkPlotPC *pc, 
						 GtkPSFont *psfont,
						 gint height);
static void pssetdash				(GtkPlotPC *pc, 
						 gdouble offset,
						 gdouble *values,
						 gint num_values);

static void ps_reencode_font			(FILE *file, char *fontname);

/*********************************************************************/
static GtkPlotPCClass *parent_class = NULL;

GtkType
gtk_plot_ps_get_type (void)
{
  static GtkType pc_type = 0;

  if (!pc_type)
    {
      GtkTypeInfo pc_info =
      {
        "GtkPlotPS",
        sizeof (GtkPlotPS),
        sizeof (GtkPlotPSClass),
        (GtkClassInitFunc) gtk_plot_ps_class_init,
        (GtkObjectInitFunc) gtk_plot_ps_init,
        /* reserved 1*/ NULL,
        /* reserved 2 */ NULL,
        (GtkClassInitFunc) NULL,
      };

      pc_type = gtk_type_unique (GTK_TYPE_PLOT_PC, &pc_info);
    }
  return pc_type;
}

static void
gtk_plot_ps_init (GtkPlotPS *ps)
{
  ps->psname = NULL;
  ps->gsaved = FALSE;
}


static void
gtk_plot_ps_class_init (GtkPlotPSClass *klass)
{
  GtkObjectClass *object_class;
  GObjectClass *gobject_class;
  GtkPlotPCClass *pc_class;

  parent_class = gtk_type_class (gtk_plot_pc_get_type ());

  object_class = (GtkObjectClass *) klass;
  gobject_class = (GObjectClass *) klass;
  pc_class = (GtkPlotPCClass *) klass;

  pc_class->init = psinit;
  pc_class->leave = psleave;
  pc_class->set_viewport = pssetviewport;
  pc_class->gsave = psgsave;
  pc_class->grestore = psgrestore;
  pc_class->set_color = pssetcolor;
  pc_class->set_dash = pssetdash;
  pc_class->set_lineattr = pssetlineattr;
  pc_class->draw_point = psdrawpoint;
  pc_class->draw_line = psdrawline;
  pc_class->draw_lines = psdrawlines;
  pc_class->draw_rectangle = psdrawrectangle;
  pc_class->draw_polygon = psdrawpolygon;
  pc_class->draw_circle = psdrawcircle;
  pc_class->draw_ellipse = psdrawellipse;
  pc_class->set_font = pssetfont;
  pc_class->draw_string = psdrawstring;

  object_class->destroy = gtk_plot_ps_destroy;
}

static void
gtk_plot_ps_destroy(GtkObject *object)
{
  GtkPlotPS *ps;

  ps = GTK_PLOT_PS(object);

  if(ps->psname){
    g_free(ps->psname);
    ps->psname = NULL;
  }
}

static void
gtk_plot_ps_construct (GtkPlotPS *ps,
		       const gchar *psname,
		       gint orientation,
		       gint epsflag,
		       gint page_size,
		       gdouble scalex,
		       gdouble scaley)
{
    gint width, height;

    ps->psname = g_strdup(psname);
    ps->orientation = orientation;
    ps->epsflag = epsflag;
    ps->page_size = page_size;
    ps->scalex = scalex;
    ps->scaley = scaley;

    switch (page_size){
    case GTK_PLOT_LEGAL:
        width = GTK_PLOT_LEGAL_W;
        height = GTK_PLOT_LEGAL_H;
        break;
    case GTK_PLOT_A4:
        width = GTK_PLOT_A4_W;
        height = GTK_PLOT_A4_H;
        break;
    case GTK_PLOT_EXECUTIVE:
        width = GTK_PLOT_EXECUTIVE_W;
        height = GTK_PLOT_EXECUTIVE_H;
        break;
    case GTK_PLOT_LETTER:
    default:
        width = GTK_PLOT_LETTER_W;
        height = GTK_PLOT_LETTER_H;
    }

    gtk_plot_ps_set_size(ps, GTK_PLOT_PSPOINTS, width, height);
}

GtkObject *
gtk_plot_ps_new (const gchar *psname,
		 gint orientation,
		 gint epsflag,
		 gint page_size,
		 gdouble scalex,
		 gdouble scaley)
{
    GtkObject *object;
    GtkPlotPS *ps;

    object = gtk_type_new(gtk_plot_ps_get_type());

    ps = GTK_PLOT_PS(object);

    gtk_plot_ps_construct(ps, psname, orientation, epsflag, page_size, scalex, scaley);

    return (object);
}

GtkObject *
gtk_plot_ps_new_with_size (const gchar *psname,
			   gint orientation,
			   gint epsflag,
			   gint units,
			   gdouble width, gdouble height,
			   gdouble scalex, gdouble scaley)
{
    GtkObject *object;
    GtkPlotPS *ps;

    object = gtk_type_new(gtk_plot_ps_get_type());

    ps = GTK_PLOT_PS(object);

    gtk_plot_ps_construct_with_size (ps, psname, orientation, epsflag, units, 
				     width, height, scalex, scaley);

    return object;
}

void
gtk_plot_ps_construct_with_size (GtkPlotPS *ps,
				 const gchar *psname,
				 gint orientation,
				 gint epsflag,
				 gint units,
				 gdouble width, gdouble height,
				 gdouble scalex, gdouble scaley)
{
    gtk_plot_ps_construct(ps, psname, orientation, epsflag, GTK_PLOT_CUSTOM, scalex, scaley);
    gtk_plot_ps_set_size(ps, units, width, height);
}

void
gtk_plot_ps_set_size (GtkPlotPS *ps,
		      gint units,
		      gdouble width,
		      gdouble height)
{
    ps->units = units;
    ps->width = width;
    ps->height = height;

    switch(units){
    case GTK_PLOT_MM:
        ps->page_width = (gdouble)width * 2.835;
        ps->page_height = (gdouble)height * 2.835;
        break;
    case GTK_PLOT_CM:
        ps->page_width = width * 28.35;
        ps->page_height = height * 28.35;
        break;
    case GTK_PLOT_INCHES:
        ps->page_width = width * 72;
        ps->page_height = height * 72;
        break;
    case GTK_PLOT_PSPOINTS:
    default:
        ps->page_width = width;
        ps->page_height = height;
    }

    if (ps->orientation == GTK_PLOT_PORTRAIT)
	gtk_plot_pc_set_viewport(GTK_PLOT_PC(ps), ps->page_width, ps->page_height);
    else
	gtk_plot_pc_set_viewport(GTK_PLOT_PC(ps), ps->page_height, ps->page_width);
}

void
gtk_plot_ps_set_scale (GtkPlotPS *ps,
		       gdouble scalex,
		       gdouble scaley)
{
    ps->scalex = scalex;
    ps->scaley = scaley; 
}

static void pssetviewport (GtkPlotPC *pc,
			   gdouble w, gdouble h) 
{

}

static void pssetlineattr (GtkPlotPC *pc, 
			   gfloat line_width,
			   GdkLineStyle line_style,
			   GdkCapStyle cap_style,
			   GdkJoinStyle join_style)
{
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    fprintf(psout, "%g slw\n", line_width);
    fprintf(psout, "%d slc\n", abs(cap_style - 1));
    fprintf(psout, "%d slj\n", join_style);

    if (line_style == 0)
	fprintf(psout,"[] 0 sd\n");  /* solid line */
}

static void 
pssetdash (GtkPlotPC *pc,
	   gdouble offset, 
	   gdouble *values,
	   gint num_values)
{
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    switch(num_values){
    case 0:
        fprintf(psout,"[] 0 sd\n");
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

void psleave (GtkPlotPC *pc)
{
    fprintf(GTK_PLOT_PS(pc)->psfile, "showpage\n");
    fprintf(GTK_PLOT_PS(pc)->psfile, "%%%%Trailer\n");
    fprintf(GTK_PLOT_PS(pc)->psfile, "%%%%EOF\n");
    fclose(GTK_PLOT_PS(pc)->psfile);
    gretl_pop_c_numeric_locale();
}

gboolean psinit (GtkPlotPC *pc)
{
    time_t now;
    FILE *psout;
    GtkPlotPS *ps;

    now = time(NULL);
    
    gretl_push_c_numeric_locale();

    ps = GTK_PLOT_PS(pc);
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

    if (ps->orientation == GTK_PLOT_PORTRAIT)
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
    /* Don't reencode the Symbol font, as it doesn't work in latin1 encoding.
     * Instead, just define Symbol-latin1 to be the same as Symbol. */
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

static void pssetcolor(GtkPlotPC *pc, const GdkColor *color)
{
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    fprintf(psout, "%g %g %g setrgbcolor\n",
	    (gdouble) color->red / 65535.0,
	    (gdouble) color->green / 65535.0,
	    (gdouble) color->blue / 65535.0);
}

static void
psdrawpoint(GtkPlotPC *pc, gdouble x, gdouble y)
{
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    fprintf(psout, "n\n");
    fprintf(psout, "%g %g m\n", x, y);
    fprintf(psout, "%g %g l\n", x, y);
    fprintf(psout, "s\n");
}

static void
psdrawlines(GtkPlotPC *pc, GtkPlotPoint *points, gint numpoints)
{
    gint i;
    FILE *psout = GTK_PLOT_PS(pc)->psfile;
 
    fprintf(psout,"n\n");
    fprintf(psout,"%g %g m\n", points[0].x, points[0].y);
    for(i = 1; i < numpoints; i++)
        fprintf(psout,"%g %g l\n", points[i].x, points[i].y);

    fprintf(psout,"s\n");
}

static void
psdrawpolygon(GtkPlotPC *pc, gboolean filled, GtkPlotPoint *points, gint numpoints)
{
    gint i;
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    fprintf(psout,"n\n");
    fprintf(psout,"%g %g m\n", points[0].x, points[0].y);
    for(i = 1; i < numpoints; i++)
	fprintf(psout,"%g %g l\n", points[i].x, points[i].y);

    if(filled)
	fprintf(psout,"f\n");
    else
	fprintf(psout,"cp\n");

    fprintf(psout,"s\n");
}

static void psdrawline(GtkPlotPC *pc, gdouble x0, gdouble y0, gdouble xf, gdouble yf)
{
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    fprintf(psout, "%g %g m\n", x0, y0);
    fprintf(psout, "%g %g l\n", xf, yf);
    fprintf(psout, "s\n");
}

static void
psdrawrectangle(GtkPlotPC *pc, gboolean filled, 
                gdouble x, gdouble y, gdouble width, gdouble height)
{
    GtkPlotPoint point[4];

    point[0].x = x;
    point[0].y = y;
    point[1].x = x + width;
    point[1].y = y;
    point[2].x = x + width;
    point[2].y = y + height;
    point[3].x = x;
    point[3].y = y + height;

    psdrawpolygon(pc, filled, point, 4);
}

static void
psdrawcircle(GtkPlotPC *pc, gboolean filled, gdouble x, gdouble y, gdouble size)
{
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    fprintf(psout,"n %g %g %g %g 0 360 ellipse\n", 
	    x, y, size / 2., size / 2.);

    if(filled)
	fprintf(psout,"f\n");

    fprintf(psout,"s\n");
}

static void
psdrawellipse(GtkPlotPC *pc, 
              gboolean filled, 
              gdouble x, gdouble y, 
              gdouble width, gdouble height)
{
    FILE *psout = GTK_PLOT_PS(pc)->psfile;

    fprintf(psout,"n %g %g %g %g 0 360 ellipse\n", 
	    x+width/2., y+height/2., 
	    width/2., height/2.);

    if(filled)
	fprintf(psout,"f\n");

    fprintf(psout,"s\n");
}

static void
psoutputstring (GtkPlotPC *pc,
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
  FILE *out = GTK_PLOT_PS(pc)->psfile;
  gchar *mbs, *c;
  

  p = wstring;
  
  if (psfont->i18n_latinfamily) {
    while (*p) {
      code = (*p <= 0x7f) ? 1 : 2;
      if (curcode && curcode != code)
	fprintf(out, "%c %s\n", end[curcode], addstring);
      if (curcode != code) {
	pssetfont(pc, fonts[code], height);
	fputc(begin[code], out);
      }
      
      curcode = code;

      wcs[0] = *p++;
      wcs[1] = 0;
      c = mbs = gdk_wcstombs(wcs);
      
      if (code == 2) {
	while (*c)
	  fprintf(out, "%02x", (unsigned char)(*c++));
      } else {
	if (*c == '(' || *c == ')')
	  fputc('\\', out);
	fputc(*c, out);
      }
      
      g_free(mbs);
    }
  } else {
    c = mbs = gdk_wcstombs(wstring);
    
    while (*c) {
      if (curcode == 0) {
	pssetfont(pc, psfont, height);
	fputc(begin[1], out);
	curcode = 1;
      }

      if (*c == '(' || *c == ')')
	fputc('\\', out);
      fputc(*c++, out);
    }

    g_free(mbs);
  }

  if (curcode)
    fprintf(out, "%c %s\n", end[curcode], addstring);
}

static void
psdrawstring	(GtkPlotPC *pc,
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

  if (text == NULL || strlen(text) == 0) return;

  psout = GTK_PLOT_PS(pc)->psfile;

  gtk_psfont_get_families(&family, &numf);
  base_psfont = psfont = gtk_psfont_get_by_name(font);
  italic = psfont->italic;
  bold = psfont->bold;

  currfont = psfont->psname;

  if (psfont->i18n_latinfamily) {
    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily, italic,
					    bold);
  }

  gtk_plot_text_get_area(text, angle, justification, font, font_height,
                         &tx, &ty, &width, &height);

  tx += x;
  ty += y;
  if(!transparent){
    pssetcolor(pc, bg);
    gtk_plot_pc_draw_rectangle(pc,
                         TRUE,
                         tx - border_space, ty - border_space,
                         width + 2*border_space, height + 2*border_space);
  }
/* border */

  pssetcolor(pc, fg);
  pssetdash(pc, 0, NULL, 0);
  pssetlineattr(pc, border_width, 0, 0, 0);
 
  switch(border){
    case GTK_PLOT_BORDER_SHADOW:
      psdrawrectangle(pc,
                         TRUE, 
                         tx - border_space + shadow_width,
                         ty + height + border_space, 
                         width + 2 * border_space, shadow_width);
      psdrawrectangle(pc,
                         TRUE, 
                         tx + width + border_space, 
                         ty - border_space + shadow_width, 
                         shadow_width, height + 2 * border_space);
    case GTK_PLOT_BORDER_LINE: 
      psdrawrectangle(pc,
                         FALSE, 
                         tx - border_space, ty - border_space, 
                         width + 2*border_space, height + 2*border_space);
    case GTK_PLOT_BORDER_NONE:
    default:
        break;
  }


  gtk_plot_text_get_size(text, angle, psfont->psname, font_height, 
                         &twidth, &theight, &tascent, &tdescent);

  if(angle == 90 || angle == 270) angle = 360 - angle;
  psgsave(pc);
  fprintf(psout, "%d %d translate\n", x, y);
  fprintf(psout, "%d rotate\n", angle);

  fprintf(psout, "0 0 m\n");
  fprintf(psout, "1 -1 sc\n");

  
  if (psfont->i18n_latinfamily)
    special = TRUE;

  c = text;
  while(c && *c != '\0' && *c != '\n') {
     if(*c == '\\'){
         c++;
         switch(*c){
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

  if(special){
    switch (justification) {
      case GTK_JUSTIFY_LEFT:
        break;
      case GTK_JUSTIFY_RIGHT:
        if(angle == 0 || angle == 180)
               fprintf(psout, "%d JR\n", twidth);
        else
               fprintf(psout, "%d JR\n", theight);
        break;
      case GTK_JUSTIFY_CENTER:
      default:
        if(angle == 0 || angle == 180)
               fprintf(psout, "%d JC\n", twidth);
        else
               fprintf(psout, "%d JC\n", theight);
        break;
    }
  } else {
    pssetfont(pc, psfont, font_height);
    
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

    psgrestore(pc);  
    fprintf(psout, "n\n");  
    return;
  }

  i = strlen(text) + 2;
  curstr = g_malloc0(sizeof(GdkWChar) * i);
  aux = wtext = g_malloc0(sizeof(GdkWChar) * i);
  gdk_mbstowcs(wtext, text, i - 1);

  scale = font_height;
  curcnt = 0;

  while(aux && *aux != '\0' && *aux != '\n') {
     if(*aux == '\\'){
         aux++;
         switch(*aux){
           case '0': case '1': case '2': case '3':
           case '4': case '5': case '6': case '7': case '9':
                  curstr[curcnt] = 0;
		  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
				 curstr, "show");
                  curcnt = 0;
                  psfont = gtk_psfont_get_by_family((gchar *)g_list_nth_data(family, *aux-'0'), italic, bold);
		  aux++;
                  break;
           case '8':case 'g':
                  curstr[curcnt] = 0;
                  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
				 curstr, "show");
                  curcnt = 0;
                  psfont = gtk_psfont_get_by_family("Symbol", italic, bold);
                  aux++;
                  break;
           case 'B':
                  curstr[curcnt] = 0;
                  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
				 curstr, "show");
                  curcnt = 0;
  		  bold = TRUE;
                  psfont = gtk_psfont_get_by_family(psfont->family, italic, bold);
		  if (psfont->i18n_latinfamily)
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily, italic, bold);
		  aux++;
                  break;
           case 'x':
                  xaux = aux + 1;
                  for (i=0; i<3; i++){
		    if (xaux[i] >= '0' && xaux[i] <= '9')
		      num[i] = xaux[i];
		    else
		      break;
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
                  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
				 curstr, "show");
                  curcnt = 0;
		  italic = TRUE;
                  psfont = gtk_psfont_get_by_family(psfont->family, italic, bold);
		  if (psfont->i18n_latinfamily)
		    latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily, italic, bold);
		  aux++;
                  break;
           case 's':case '_':
                  curstr[curcnt] = 0;
                  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
				 curstr, "show");
                  curcnt = 0;
                  scale = 0.6 * font_height;
		  offset -= (gint)scale / 2;
                  fprintf(psout, "0 %d rmoveto\n", -((gint)scale / 2));
                  aux++;
                  break;
           case 'S':case '^':
                  curstr[curcnt] = 0;
                  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
				 curstr, "show");
                  curcnt = 0;
                  scale = 0.6 * font_height;
		  offset += 0.5*font_height;
                  fprintf(psout, "0 %d rmoveto\n", (gint)(0.5*font_height));
                  aux++;
                  break;
           case 'N':
                  curstr[curcnt] = 0;
                  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
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
                  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
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
		  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
				 bkspchar,
				 "stringwidth pop 0 exch neg exch rmoveto");
                  aux++;
                  break;
           case '-':
                  curstr[curcnt] = 0;
		  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
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
		  psoutputstring(pc, psfont, latin_psfont, (gint)scale,
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
     } else {
       if(aux && *aux != '\0' && *aux != '\n'){
                curstr[curcnt++] = *aux;
		lastchar = aux;
                aux++;

       }
     }
  }
  curstr[curcnt] = 0;
  psoutputstring(pc, psfont, latin_psfont, (gint)scale, curstr, "show");

  psgrestore(pc);  
  fprintf(psout, "n\n");  


  g_free(wtext);
  g_free(curstr);
}

static void
pssetfont(GtkPlotPC *pc, GtkPSFont *psfont, gint height)
{
  FILE *psout = GTK_PLOT_PS(pc)->psfile;

  if (psfont->i18n_latinfamily && psfont->vertical)
    fprintf(psout,
	    "/%s ff [0 1 -1 0 0 0.3] makefont [%d 0 0 %d 0 0] makefont sf\n",
	    psfont->psname, height, height);
  else
    fprintf(psout, "/%s-latin1 ff %g scf sf\n", psfont->psname, (double)height);

}


static void
psgsave(GtkPlotPC *pc)
{
  GtkPlotPS *ps;
  FILE *psout;

  ps = GTK_PLOT_PS(pc);

  psout = ps->psfile;

  fprintf(psout,"gsave\n");
  ps->gsaved = TRUE;
}

static void
psgrestore(GtkPlotPC *pc)
{
  GtkPlotPS *ps;
  FILE *psout;

  ps = GTK_PLOT_PS(pc);

  psout = ps->psfile;

  if(!ps->gsaved) return;

  fprintf(psout,"grestore\n");
  ps->gsaved = FALSE;
}

static void
gtk_plot_text_get_size(const gchar *text, gint angle, 
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
  aux = wtext = g_malloc0(sizeof(GdkWChar) * i);
  gdk_mbstowcs(wtext, text, i - 1);

  while(aux && *aux != '\0' && *aux != '\n'){
   if(*aux == '\\'){
     aux++;
     switch(*aux){
       case '0': case '1': case '2': case '3':
       case '4': case '5': case '6': case '7': case '9':
           psfont = gtk_psfont_get_by_family((gchar *)g_list_nth_data(family, *aux-'0'), italic, bold);
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
           if(latin_font){
             gdk_font_unref(latin_font);
             latin_font = NULL;
           }
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
           if(latin_font){
             gdk_font_unref(latin_font);
             latin_font = NULL;
           }
           if (psfont->i18n_latinfamily) {
             latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily,
                                                      italic, bold);
             latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
           }
           aux++;
           break;
       case 'S': case '^':
           fontsize = (int)((gdouble)fontsize * 0.6 + 0.5);
           gdk_font_unref(font);
           font = gtk_psfont_get_gdkfont(psfont, fontsize);
           if(latin_font){
             gdk_font_unref(latin_font);
             latin_font = NULL;
           }
           if (psfont->i18n_latinfamily)
             latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);

	   y -= font->ascent;
           aux++;
           break;
       case 's': case '_':
           fontsize = (int)((gdouble)fontsize * 0.6 + 0.5);
           gdk_font_unref(font);
           font = gtk_psfont_get_gdkfont(psfont, fontsize);
           if(latin_font){
             gdk_font_unref(latin_font);
             latin_font = NULL;
           }
           if (psfont->i18n_latinfamily)
             latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);

           y += font->descent;
           aux++;
           break;
       case '+':
           fontsize += 3;
           gdk_font_unref(font);
           font = gtk_psfont_get_gdkfont(psfont, fontsize);
           if(latin_font){
             gdk_font_unref(latin_font);
             latin_font = NULL;
           }
           if (psfont->i18n_latinfamily)
             latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);

           aux++;
           break;
       case '-':
           fontsize -= 3;
           gdk_font_unref(font);
           font = gtk_psfont_get_gdkfont(psfont, fontsize);
           if(latin_font){
             gdk_font_unref(latin_font);
             latin_font = NULL;
           }
           if (psfont->i18n_latinfamily)
             latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);

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
           if(latin_font){
             gdk_font_unref(latin_font);
             latin_font = NULL;
           }
           if (psfont->i18n_latinfamily) {
             latin_psfont = gtk_psfont_get_by_family(psfont->i18n_latinfamily,
                                                      italic, bold);
             latin_font = gtk_psfont_get_gdkfont(latin_psfont, fontsize);
           }

           aux++;
           break;
       case 'b':
	   if(lastchar){
	     gtk_psfont_get_char_size(psfont, font, latin_font, *lastchar, &w,
				      NULL, NULL);
	     x -= w;
	     
	     if (lastchar == wtext)
	       lastchar = NULL;
	     else
	       lastchar--;
	   } else {
	     gtk_psfont_get_char_size(psfont, font, latin_font, 'X', &w, NULL, NULL);
	     x -= w;
           }
           aux++;
	   break;
       default:
           if(aux && *aux != '\0' && *aux !='\n'){
	     gtk_psfont_get_char_size(psfont, font, latin_font, *aux, &w, &a, &d);
	     x += w;
	     lastchar = aux;
             aux++;
           }
           break;
     }
   } else {
     if(aux && *aux != '\0' && *aux != '\n'){
       gtk_psfont_get_char_size(psfont, font, latin_font, *aux, &w, &a, &d);
       x += w;
       lastchar = aux;
       aux++;
       if(x > old_width) old_width = x;
       if(y + d - y0 > *descent) *descent = y + d - y0;
       if(y0 - y + a > *ascent) *ascent = y0 - y + a;
     }
   }
  }

  old_height = *ascent + *descent;
  *width = old_width;
  *height = old_height;
  if(angle == 90 || angle == 270)
    {
      *width = old_height;
      *height = old_width;
    }

  g_free(wtext);
  gdk_font_unref(font);
}

static void
gtk_plot_text_get_area(const gchar *text, gint angle, GtkJustification just, 
                       const gchar *font, gint font_height,
                       gint *x, gint *y,
                       gint *width, gint *height)
{
  gint ascent, descent;

  if(text == NULL) return;

  gtk_plot_text_get_size(text, angle, font, 
                         font_height, 
                         width, height, &ascent, &descent);

  *x = 0;
  *y = 0;

  switch(just){
    case GTK_JUSTIFY_LEFT:
      switch(angle){
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
      switch(angle){
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
      switch(angle){
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








