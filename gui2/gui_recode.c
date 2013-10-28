/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "gretl.h"
#include <glib/gstdio.h>

static int seven_bit_string (const unsigned char *s)
{
    while (*s) {
	if (*s > 127) return 0;
	s++;
    }
    return 1;
}

/* See if @fname works as an argument to g_fopen in read mode. If so,
   return 0 and ignore the @fconv argument. If not, try re-encoding
   into @fconv and if that works then return 0. If we can't get
   g_fopen to work either way, display an error message and return
   E_FOPEN.
*/

int validate_filename_for_glib (const gchar *fname, gchar **fconv)
{
    FILE *fp = g_fopen(fname, "r");
    int err = 0;

    if (fp == NULL) {
	gsize bytes;

	*fconv = NULL;

#ifdef G_OS_WIN32
	if (!g_utf8_validate(fname, -1, NULL)) {
	    /* The Glib filename encoding is UTF-8 on Windows */
	    *fconv = g_locale_to_utf8(fname, -1, NULL, &bytes, NULL);
	}
#else
	if (g_utf8_validate(fname, -1, NULL)) {
	    /* *nix: maybe we shouldn't be using UTF-8? */
	    *fconv = g_filename_from_utf8(fname, -1, NULL, &bytes, NULL);
	}
#endif
	if (*fconv != NULL) {
	    fp = g_fopen(*fconv, "r");
	    if (fp == NULL) {
		g_free(*fconv);
	    }
	}
    }

    if (fp != NULL) {
	fclose(fp);
    } else {
	err = E_FOPEN;
	errbox(_("Can't open %s for reading"), fname);
    }

    return err;
}

/* This is used for converting the UTF-8 datafile name
   inside a gretl session file to the locale.
*/

gchar *my_filename_from_utf8 (char *fname)
{
    const gchar *cset;
    GError *err = NULL;
    gsize bytes;
    gchar *tmp = NULL;

    if (seven_bit_string((unsigned char *) fname)) {
	return fname;
    }

    if (g_get_charset(&cset)) {
	/* locale uses UTF-8, OK */
	return fname;
    }

#ifdef G_OS_WIN32
    /* don't use g_filename_from_utf8 because on Windows
       GLib uses UTF-8 for filenames and no conversion
       will take place */
    tmp = g_locale_from_utf8(fname, -1, NULL, &bytes, &err);
#else
    tmp = g_filename_from_utf8(fname, -1, NULL, &bytes, &err);
#endif

    if (err) {
	errbox(err->message);
	g_error_free(err);
    } else {
	strcpy(fname, tmp);
    }

    g_free(tmp);

    return fname;
}

/* returns new copy of fname, converted if need be */

gchar *my_filename_to_utf8 (const char *fname)
{
    GError *err = NULL;
    gchar *ret = NULL;

    if (g_utf8_validate(fname, -1, NULL)) {
	ret = g_strdup(fname);
    } else {
	/* On Windows, with GTK >= 2.6, the GLib filename
	   encoding is UTF-8; however, filenames coming from
	   a native Windows source will be in the
	   locale charset 
	*/
	gsize bytes;

#ifdef G_OS_WIN32
	ret = g_locale_to_utf8(fname, -1, NULL, &bytes, &err);
#else
	ret = g_filename_to_utf8(fname, -1, NULL, &bytes, &err);
#endif
    }

    if (err) {
	errbox(err->message);
	g_error_free(err);
    } 

    return ret;
}

gchar *my_locale_from_utf8 (const gchar *src)
{
    const gchar *cset;
    gsize bytes;
    GError *err = NULL;
    gchar *ret = NULL;

    if (g_get_charset(&cset)) {
	/* g_get_charset returns TRUE if the returned 
	   charset is UTF-8 */ 
	return g_strdup(src);
    } 

    ret = g_locale_from_utf8(src, -1, NULL, &bytes, &err);

    if (err) {
	errbox(err->message);
	g_error_free(err);
    }

    return ret;
}

/* Used when, e.g. loading a script into a GTK window: the
   script might be in a locale encoding other than UTF-8
   (if it was created in a third-party editor).
*/

static gchar *real_my_locale_to_utf8 (const gchar *src,
				      int starting)
{
    static int errcount;
    const gchar *charset = NULL;
    gsize bytes;
    GError *err = NULL;
    gchar *ret;

    if (g_get_charset(&charset)) {
	/* in a UTF-8 locale */
	if (starting) {
	    errcount = 0;
	    ret = g_convert(src, -1, "UTF-8", "ISO-8859-15", 
			    NULL, &bytes, &err);
	} else if (errcount == 0) {
	    ret = g_convert(src, -1, "UTF-8", "ISO-8859-15", 
			    NULL, &bytes, &err);
	} else {
	    ret = g_convert(src, -1, "UTF-8", "ISO-8859-15", 
			    NULL, &bytes, NULL);
	}	
    } else {
	if (starting) {
	    errcount = 0;
	    ret = g_locale_to_utf8(src, -1, NULL, &bytes, &err);
	} else if (errcount == 0) {
	    ret = g_locale_to_utf8(src, -1, NULL, &bytes, &err);
	} else {
	    ret = g_locale_to_utf8(src, -1, NULL, &bytes, NULL);
	}
    }

    if (err) {
	errbox(err->message);
	g_error_free(err);
	errcount++;
    }

    return ret;
}

/* wrappers to avoid repeated error messages where one
   will make the point
*/

gchar *my_locale_to_utf8 (const gchar *src)
{
    return real_my_locale_to_utf8(src, 1);
}

gchar *my_locale_to_utf8_next (const gchar *src)
{
    return real_my_locale_to_utf8(src, 0);
}

/* We're not totally sure what charset we're supposed to be converting
   from, but we'll make a guess!
*/

static gchar *gp_locale_to_utf8 (const gchar *src, int starting)
{
    static int errcount;
    static const char *from;
    gsize read, wrote;
    GError *err = NULL;
    gchar *ret;

    if (from == NULL) {
	if (iso_latin_version() == 2) {
#ifdef G_OS_WIN32
	    from = "CP1250";
#else
	    from = "ISO-8859-2";
#endif
	} else {
	    from = "ISO-8859-1";
	}
    } 

    if (starting) {
	errcount = 0;
	ret = g_convert(src, -1, "UTF-8", from,
			&read, &wrote, &err);
    } else if (errcount == 0) {
	ret = g_convert(src, -1, "UTF-8", from,
			&read, &wrote, &err);
    } else {
	ret = g_convert(src, -1, "UTF-8", from,
			&read, &wrote, NULL);
    }

    if (err != NULL) {
	errbox(err->message);
	g_error_free(err);
	errcount++;
    }

    return ret;
}

/* gp_locale_from_utf8: used when taking gnuplot commands
   from a GTK editor window and sending them to gnuplot
   or saving to "user file", on non-UTF-8 platforms.
*/

gchar *gp_locale_from_utf8 (const gchar *src)
{
    gsize read, wrote;
    GError *err = NULL;
    gchar *ret;

    if (gretl_is_ascii(src)) {
	return NULL;
    }

    /* let glib figure out what the target locale is */

    ret = g_locale_from_utf8(src, -1, &read, &wrote, &err);

    if (err) {
	errbox(err->message);
	g_error_free(err);
    }

    return ret;
}

static int maybe_adjust_cairo (char *line)
{
    char *s = line + 12;
    int ret = 0;

     if (!strncmp(s, "cairo", 5)) {
	if (gnuplot_png_terminal() != GP_PNG_CAIRO) {
	    /* drop back to non-cairo PNG term */
	    shift_string_left(s + 5, 5);
	    ret = 1;
	}
    } else if (gnuplot_png_terminal() == GP_PNG_CAIRO &&
	       strlen(line) < 512 - 5) {
	/* substitute the preferred cairo PNG term */
	char *p, tmp[512];

	strcpy(tmp, "set term pngcairo");
	s += strspn(s, " ");
	if (!strncmp(s, "truecolor ", 10)) {
	    /* invalid */
	    s += 10;
	}
	if (!strncmp(s, "font", 4)) {
	    /* attempt to fix up old gnuplot font spec */
	    char fname[32], fsize[6];

	    s += 4;
	    s += strspn(s, " ");
	    if (*s != '"') {
		if (sscanf(s, "%31s %5s", fname, fsize) == 2) {
		    strcat(tmp, " font \"");
		    strcat(tmp, fname);
		    strcat(tmp, ",");
		    strcat(tmp, fsize);
		    strcat(tmp, "\" ");
		}
	    }
	    p = strstr(s, " size ");
	    if (p != NULL) {
		strcat(tmp, p + 1);
	    } else {
		strncat(tmp, "\n", 1);
	    }
	} else {
	    /* hope for the best */
	    strcat(tmp, s);
	}
	strcpy(line, tmp);
	ret = 1;
    }

    return ret;
}

static void do_fix_xrange (char *line, int fix, FILE *fp)
{
    double xmin, xmax;
    int n;

    gretl_push_c_numeric_locale();

    n = sscanf(line, "set xrange [%lf:%lf]", &xmin, &xmax);

    /* note: if we fail here, we'll suppress the printing of
       the xrange specification */

    if (n == 2) {
	double offset = 946684800;

	if (fix == 1) {
	    /* old gnuplot to new */
	    xmin += offset;
	    xmax += offset;
	} else {
	    /* new gnuplot to old */
	    xmin -= offset;
	    xmax -= offset;
	}
	fprintf(fp, "set xrange [%.12g:%.12g]\n", xmin, xmax);
    }

    gretl_pop_c_numeric_locale();
}

/* Backward compatibility for gnuplot command files as saved in
   sessions: if the file is non-ascii and non-UTF-8, convert to UTF-8,
   since we have now (2008-01) standardized on UTF-8 as the encoding
   for all gnuplot files that are "internal" to gretl.

   While we're at it, check for a couple of other things.  If we have
   a "set missing" or "set datafile missing" line, ensure that it's in
   sync with the installed gnuplot.  Also, do we have a commented-out
   "set term" line?  If so, this may be out of date with respect to
   current gnuplot, so we'll strip it out.
*/

int maybe_rewrite_gp_file (const char *fname)
{
    FILE *fin, *fout;
    gchar *trbuf, *modname = NULL;
    char line[512];
    double gpver;
    int modified = 0;
    int recoded = 0;
    int fix_xrange = 0;
    int err = 0;

    fin = gretl_fopen(fname, "r");
    if (fin == NULL) {
	return 1;
    }

    modname = g_strdup_printf("%s.tr", fname);

    fout = gretl_fopen(modname, "w");
    if (fout == NULL) {
	fclose(fin);
	g_free(modname);
	return 1;
    }

    gpver = gnuplot_version();

    while (fgets(line, sizeof line, fin)) {
	int modline = 0;
	
	if (!strncmp(line, "set missing", 11)) {
	    fputs("set datafile missing \"?\"\n", fout);
	    modline = 1;
	} else if (!strncmp(line, "# set term", 10)) {
	    /* skip it */
	    modline = 1;
	} else if (!strncmp(line, "set term png", 12)) {
	    if (maybe_adjust_cairo(line)) {
		fputs(line, fout);
		modline = 1;
	    }
	} else if (!strncmp(line, "set xdata time", 14)) {
	    int zy2000 = strstr(line, "ZERO_YEAR=2000") != NULL;

	    if (gpver >= 4.7 && zy2000) {
		fputs("set xdata time\n", fout);
		modline = 1;
		fix_xrange = 1;
	    } else if (gpver < 4.7 && !zy2000) {
		fputs("set xdata time # ZERO_YEAR=2000\n", fout);
		modline = 1;
		fix_xrange = 2;
	    }
	} else if (!strncmp(line, "set xrange", 10) && fix_xrange) {
	    do_fix_xrange(line, fix_xrange, fout);
	    modline = 1;
	}

	if (!modline && !g_utf8_validate(line, -1, NULL)) {
	    trbuf = gp_locale_to_utf8(line, !recoded);
	    if (trbuf != NULL) {
		fputs(trbuf, fout);
		g_free(trbuf);
	    } 
	    modline = recoded = 1;
	}

	if (modline) {
	    modified = 1;
	} else {
	    /* pass old line through */
	    fputs(line, fout);
	}
    }

    fclose(fin);
    fclose(fout);

    if (modified) {
	err = copyfile(modname, fname);
    }

    gretl_remove(modname);
    g_free(modname);

    return err;
}

