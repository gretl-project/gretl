/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

/*  htmlprint.h for gretl */ 

#ifndef HTMLPRINT_H
#define HTMLPRINT_H

#ifdef USE_GTKHTML
# include <gtkhtml/gtkhtml.h>
#else
 typedef void GtkHTML;
 typedef void GtkHTMLStream;
#define GTK_IS_HTML(w)  FALSE
#endif

typedef struct { 
    GtkHTML *w;              /* widget */
    GtkHTMLStream *strm;     /* stream */
    FILE *fp;                /* for printing html to file */
} html_t;

#ifdef USE_GTKHTML
int h_bufopen (html_t *htm);
#endif
int h_fopen (html_t *htm, char *fname);
void h_bufclose (html_t *htm);
void h_printmodel (const MODEL *pmod, 
		   const DATAINFO *pdinfo, 
		   html_t *htm);

#endif /* HTMLPRINT_H */





