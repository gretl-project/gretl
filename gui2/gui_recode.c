/*
 *  Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "gretl.h"

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
    gchar *trfname = NULL;
    GError *err = NULL;
    gsize bytes;

    if (seven_bit_string((unsigned char *) fname)) {
	return fname;
    }

    trfname = g_filename_from_utf8(fname, -1, NULL, &bytes, &err);

    if (err != NULL) {
	errbox(err->message);
	g_error_free(err);
    } else {
	strcpy(fname, trfname);
    }

    g_free(trfname);

    return fname;
}

#define TR_DEBUG 0

static gchar *
real_locale_from_utf8 (const gchar *src, int force)
{
    gchar *trstr;
    gsize bytes;
    GError *err = NULL;
    const gchar *cset = NULL;

#if TR_DEBUG
    gchar *msg;
    int u = g_get_charset(&cset);
    msg = g_strdup_printf("real_locale_from_utf8: force=%d, "
			  "g_get_charset returned %d and gave "
			  "'%s'\n", force, u, 
			  (cset != NULL)? cset : "NULL");
    infobox(msg);
    g_free(msg);

    if (!force && u) {
	return g_strdup(src);
    }
#else
    if (!force && g_get_charset(&cset)) {
	/* According to the glib manual, g_get_charset returns TRUE if 
	   the returned charset is UTF-8 */ 
	return g_strdup(src);
    }
#endif

    trstr = g_locale_from_utf8(src, -1, NULL, &bytes, &err);

    if (err != NULL) {
	if (cset != NULL) {
	    errbox("g_locale_from_utf8 failed for charset '%s'", cset);
	} else {
	    errbox("g_locale_from_utf8 failed; so did g_get_charset");
	}
	g_error_free(err);
    }

    return trstr;
}

gchar *my_locale_from_utf8 (const gchar *src)
{
    return real_locale_from_utf8(src, 0);
}

gchar *force_locale_from_utf8 (const gchar *src)
{
    return real_locale_from_utf8(src, 1);
}

#define FNAME_DEBUG 0

gchar *my_filename_to_utf8 (char *fname)
{
    gchar *trfname = NULL;
    GError *err = NULL;
    gsize bytes;

#if FNAME_DEBUG
    fprintf(stderr, "my_filename_to_utf8: fname='%s'\n", fname);
    fflush(stderr);
#endif

    if (g_utf8_validate(fname, -1, NULL)) {
#if FNAME_DEBUG
	fprintf(stderr, " validates as utf8, returning fname\n");
	fflush(stderr);
#endif
	return fname;
    }

    trfname = g_filename_to_utf8(fname, -1, NULL, &bytes, &err);

    if (err != NULL) {
	errbox("g_filename_to_utf8 failed");
	g_error_free(err);
    } else {
	strcpy(fname, trfname);
#if FNAME_DEBUG
	fprintf(stderr, " converted fname='%s'\n", fname);
	fflush(stderr);
#endif
    }

    g_free(trfname);

    return fname;
}

static gchar *real_my_locale_to_utf8 (const gchar *src,
				      int starting)
{
    static int errcount;
    gchar *trstr;
    gsize bytes;
    GError *err = NULL;

    if (starting) {
	errcount = 0;
	trstr = g_locale_to_utf8(src, -1, NULL, &bytes, &err);
    } else if (errcount == 0) {
	trstr = g_locale_to_utf8(src, -1, NULL, &bytes, &err);
    } else {
	trstr = g_locale_to_utf8(src, -1, NULL, &bytes, NULL);
    }

    if (err != NULL) {
	const gchar *cset = NULL;

	g_get_charset(&cset);
	if (cset != NULL) {
	    errbox("g_locale_to_utf8 failed for charset '%s'", cset);
	} else {
	    errbox("g_locale_to_utf8 failed; so did g_get_charset");
	}
	g_error_free(err);
	errcount++;
    }

    return trstr;
}

gchar *my_locale_to_utf8 (const gchar *src)
{
    return real_my_locale_to_utf8(src, 1);
}

gchar *my_locale_to_utf8_next (const gchar *src)
{
    return real_my_locale_to_utf8(src, 0);
}

gchar *latin2_to_utf8 (const gchar *src)
{
    gsize read, wrote;

#ifdef G_OS_WIN32
    return g_convert(src, -1, "UTF-8", "CP1250",
		     &read, &wrote, NULL);
#else
    return g_convert(src, -1, "UTF-8", "ISO-8859-2",
		     &read, &wrote, NULL);
#endif
}
