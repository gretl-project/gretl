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

#ifdef ENABLE_NLS

static int seven_bit_string (const unsigned char *s)
{
    while (*s) {
	if (*s > 127) return 0;
	s++;
    }
    return 1;
}

static int seven_bit_file (const char *fname)
{
    FILE *fp;
    char line[256];
    int ascii = 1;
    
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	if (!seven_bit_string((unsigned char *) line)) {
	    ascii = 0;
	    break;
	}
    }

    fclose(fp);

    return ascii;
}

int maybe_recode_file (const char *fname)
{
    const gchar *charset;

    if (g_get_charset(&charset)) {
	/* locale uses UTF-8 */
	return 0;
    }

    if (seven_bit_file(fname)) {
	return 0;
    } else {
	FILE *fin, *fout;
	char trname[MAXLEN];
	char line[128];
	gchar *trbuf;
	int err = 0;

	fin = gretl_fopen(fname, "r");
	if (fin == NULL) {
	    return 1;
	}

	sprintf(trname, "%s.tr", fname);

	fout = gretl_fopen(trname, "w");
	if (fout == NULL) {
	    fclose(fin);
	    return 1;
	}	

	while (fgets(line, sizeof line, fin) && !err) {
	    trbuf = my_locale_from_utf8(line);
	    if (trbuf != NULL) {
		fputs(trbuf, fout);
		g_free(trbuf);
	    } else {
		err = 1;
	    }
	}

	fclose(fin);
	fclose(fout);

	if (!err) {
	    err = copyfile(trname, fname);
	    remove(trname);
	}

	return err;
    }

    return 0;
}

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
	/* is this right? */
	return fname;
    }

    tmp = g_filename_from_utf8(fname, -1, NULL, &bytes, &err);

    if (err) {
	errbox(err->message);
	g_error_free(err);
    } else {
	strcpy(fname, tmp);
    }

    g_free(tmp);

    return fname;
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

/* returns new copy of fname, converted if need be */

gchar *my_filename_to_utf8 (const char *fname)
{
    GError *err = NULL;
    gsize bytes;
    gchar *ret = NULL;

    if (g_utf8_validate(fname, -1, NULL)) {
	ret = g_strdup(fname);
    } else {
	/* On Windows, with GTK >= 2.6, the GLib filename
	   encoding is UTF-8; however, filenames coming from
	   a native Windows file dialog will be in the
	   locale charset 
	*/
#if defined(G_OS_WIN32) && GTK_MINOR_VERSION >= 6
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

static gchar *real_my_locale_to_utf8 (const gchar *src,
				      int starting)
{
    static int errcount;
    gsize bytes;
    GError *err = NULL;
    gchar *ret;

    if (starting) {
	errcount = 0;
	ret = g_locale_to_utf8(src, -1, NULL, &bytes, &err);
    } else if (errcount == 0) {
	ret = g_locale_to_utf8(src, -1, NULL, &bytes, &err);
    } else {
	ret = g_locale_to_utf8(src, -1, NULL, &bytes, NULL);
    }

    if (err) {
	errbox(err->message);
	g_error_free(err);
	errcount++;
    }

    return ret;
}

static const char *gp_cset (void)
{
    if (iso_latin_version() == 2) {
#ifdef G_OS_WIN32
	return "CP1250";
#else
	return "ISO-8859-2";
#endif
    } else {
	return "ISO-8859-1";
    } 
}

static gchar *real_gp_locale_to_utf8 (const gchar *src,
				      int starting)
{
    static int errcount;
    static const char *cset;
    gsize read, wrote;
    GError *err = NULL;
    gchar *ret;

    if (cset == NULL) {
	cset = gp_cset();
    }

    if (starting) {
	errcount = 0;
	ret = g_convert(src, -1, "UTF-8", cset,
			&read, &wrote, &err);
    } else if (errcount == 0) {
	ret = g_convert(src, -1, "UTF-8", cset,
			&read, &wrote, &err);
    } else {
	ret = g_convert(src, -1, "UTF-8", cset,
			&read, &wrote, NULL);
    }

    if (err) {
	errbox(err->message);
	g_error_free(err);
	errcount++;
    }

    return ret;
}

gchar *my_locale_to_utf8 (const gchar *src)
{
    return real_my_locale_to_utf8(src, 1);
}

gchar *my_locale_to_utf8_next (const gchar *src)
{
    return real_my_locale_to_utf8(src, 0);
}

gchar *gp_locale_to_utf8 (const gchar *src)
{
    return real_gp_locale_to_utf8(src, 1);
}

gchar *gp_locale_to_utf8_next (const gchar *src)
{
    return real_gp_locale_to_utf8(src, 0);
}

gchar *gp_locale_from_utf8 (const gchar *src)
{
    static const char *cset;
    gsize read, wrote;
    GError *err = NULL;
    gchar *ret;

    if (gretl_is_ascii(src)) {
	return NULL;
    }

    if (cset == NULL) {
	cset = gp_cset();
    }

#if 0
    fprintf(stderr, "gp_locale_from_utf8: src is not ascii, cset = '%s'\n",
	    cset);
#endif

    ret = g_convert(src, -1, cset, "UTF-8",
		    &read, &wrote, &err);

    if (err) {
	errbox(err->message);
	g_error_free(err);
    }

    return ret;
}

gchar *latin1_to_utf8 (const gchar *src)
{
    gsize read, wrote;
    GError *err = NULL;
    gchar *ret;

    ret = g_convert(src, -1, "UTF-8", "ISO-8859-1",
		    &read, &wrote, &err);

    if (err) {
	errbox(err->message);
	g_error_free(err);
    }

    return ret;
}

gchar *latin2_to_utf8 (const gchar *src)
{
    gsize read, wrote;
    GError *err = NULL;
    gchar *ret;

# ifdef G_OS_WIN32
    ret = g_convert(src, -1, "UTF-8", "CP1250",
		    &read, &wrote, &err);
# else
    ret = g_convert(src, -1, "UTF-8", "ISO-8859-2",
		    &read, &wrote, &err);
# endif

    if (err) {
	errbox(err->message);
	g_error_free(err);
    }

    return ret;
}

#endif /* ENABLE_NLS */


