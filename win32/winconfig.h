#ifndef WINCONFIG_H
#define WINCONFIG_H

/* Native language support */
#define ENABLE_NLS 1
#define PACKAGE "gretl"

/* Flag the fact that we're building a self-installer package */
#define PKGBUILD 1

/* Flag use of libcurl */
#define USE_CURL 1

/* Extra floating-point GMP routines? */
#define HAVE_MPFR 1

/* openmp? */
#define OPENMP_BUILD 1

/* Are we supporting MPI? */
#define HAVE_MPI 1

/* sse2: we'll assume this is OK */
#define USE_SSE2 1

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

/* Define if using libgsf (>= 1.14.29) */
#undef USE_GSF

/* Is LaTeX available? */
#define HAVE_LATEX 1

/* Use gtksourceview-2.0 for syntax highlighting */
#define USE_GTKSOURCEVIEW_2 1

/* Does unistd.h have getdomainname? */
#define GETDOMAINNAME 1

/* Do we have paths.h? */
/* #define HAVE_PATHS_H 1 */

/* Do we have stdint.h? */
#define HAVE_STDINT_H 1 

/* Do we have unistd.h? */
#define HAVE_UNISTD_H 1 

/* Define if the `long double' type works.  */
/* #undef HAVE_LONG_DOUBLE */

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if your processor stores words with the most significant
   byte first. */
/* #undef WORDS_BIGENDIAN */

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1

/* Audio graph support */
#define HAVE_AUDIO 1

/* Mailer (send-to) support */
#define ENABLE_MAILER 1

/* Use the GNU R shared library */
#define USE_RLIB 1

#define WINVER 0x0500 /* XP */

/* sapi speech support */
#define WIN32_SAPI 1

#endif /* WINCONFIG_H */
