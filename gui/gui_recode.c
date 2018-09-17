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
#include "gui_recode.h"

#include <glib/gstdio.h>

#ifdef G_OS_WIN32

/* Windows-specific, currently used in several gui C files:
   is it now redundant? */

gchar *my_filename_to_utf8 (const char *fname)
{
    GError *gerr = NULL;
    gsize bytes;
    gchar *ret;

    fprintf(stderr, "my_filename_to_utf8: '%s'\n", fname);

    ret = g_locale_to_utf8(fname, -1, NULL, &bytes, &gerr);

    if (gerr) {
	errbox(gerr->message);
	g_error_free(gerr);
    }

    return ret;
}

/* this variant of my_filename_to_utf8() will never
   return NULL */

gchar *filename_to_utf8_nofail (const char *fname)
{
    gchar *ret = my_filename_to_utf8(fname);

    if (ret == NULL) {
	ret = g_strdup("unknown filename");
    }

    return ret;
}

#endif

/* Used when, e.g. loading a script into a GTK window: the
   script might be in a locale encoding other than UTF-8
   if it was created in a third-party editor.
*/

static gchar *real_my_locale_to_utf8 (const gchar *src,
				      int starting)
{
    static int errcount;
    const gchar *charset = NULL;
    gsize bytes;
    GError *err = NULL;
    GError **errp = NULL;
    gchar *ret;

    if (starting) {
	errcount = 0;
    }

    /* avoid multiple repetition of error message */
    errp = errcount == 0 ? &err : NULL;

    if (g_get_charset(&charset)) {
	/* we're in a UTF-8 locale */
	ret = g_convert(src, -1, "UTF-8", "ISO-8859-15",
			NULL, &bytes, errp);
    } else {
	ret = g_locale_to_utf8(src, -1, NULL, &bytes, errp);
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

/* If need be, fix up an old, pre-cairo PNG terminal line:
   we come here only if a plot line starts "set term png"
   followed by a space rather than "cairo".
*/

static int maybe_adjust_cairo (char *line)
{
    char *s = line + 12; /* skip "set term png" */
    int ret = 0;

    if (strlen(line) < 512 - 5) {
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
		strcat(tmp, "\n");
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

    while (fgets(line, sizeof line, fin)) {
	int modline = 0;

	if (!strncmp(line, "set missing", 11)) {
	    fputs("set datafile missing \"?\"\n", fout);
	    modline = 1;
	} else if (!strncmp(line, "# set term", 10)) {
	    /* skip it */
	    modline = 1;
	} else if (!strncmp(line, "set term png ", 13)) {
	    if (maybe_adjust_cairo(line)) {
		fputs(line, fout);
		modline = 1;
	    }
	} else if (!strncmp(line, "set xdata time", 14)) {
	    if (strstr(line, "ZERO_YEAR=2000")) {
		fputs("set xdata time\n", fout);
		modline = 1;
		fix_xrange = 1;
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
