/* config.h.  Generated automatically by configure.  */
/* Configuration header file.
   Copyright (C) 1995, 1996, 1997, 1998 Free Software Foundation, Inc.

This file is part of gretl.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

#ifndef WINCONFIG_H
#define WINCONFIG_H

/* Native language support */
#define ENABLE_NLS 1
#define PACKAGE "gretl"

/* GMP library support */
#define ENABLE_GMP 1

/* X-12-ARIMA support? */
#define HAVE_X12A 1

/* TRAMO/SEATS support? */
#define HAVE_TRAMO 1

/* Define if you want GNU readline support */
#define HAVE_READLINE 1

/* Current readline? */
#define NEW_READLINE 1

/* Define if zlib is available */
#define HAVE_ZLIB 1

/* Is gnuplot available? */
#define HAVE_GNUPLOT 1

/* Use gnuplot PNG output? */
#define GNUPLOT_PNG 1

/* Is LaTeX available? */
#define HAVE_LATEX 1

/* Is Gnome available? */
/* #undef HAVE_GNOME 1 */

/* Use installed gtkextra? */
/* #undef HAVE_GTKEXTRA 1 */

/* Use gtksourceview for syntax highlighting */
#define USE_GTKSOURCEVIEW 1

/* Does unistd.h have getdomainname? */
#define GETDOMAINNAME 1

/* Do we have paths.h? */
/* #define HAVE_PATHS_H 1 */

/* Do we have unistd.h? */
#define HAVE_UNISTD_H 1 

/* Define if the `long double' type works.  */
/* #undef HAVE_LONG_DOUBLE */

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1

/* Audio graph support */
#define HAVE_AUDIO 1

/* Mailer (send-to) support */
#define ENABLE_MAILER 1

#endif /* WINCONFIG_H */
