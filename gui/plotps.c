/* functions borrowed from gtkplotps.c in gtkextra.  They are declared
   static there, but gretl wants direct access to them for the
   postscript printing of boxplots. */

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


static void ps_reencode_font (FILE *file, char *fontname);
static gboolean psinit (GtkPlotPS *ps);
static void psleave (GtkPlotPS *ps);

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
