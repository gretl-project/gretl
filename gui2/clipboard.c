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
#include "textutil.h"

#define CLIPDEBUG 0

static gchar *clipboard_buf; 

GtkTargetEntry text_targets[] = {
    { "UTF8_STRING",     0, TARGET_UTF8_STRING },
    { "STRING",          0, TARGET_STRING },
    { "TEXT",            0, TARGET_TEXT }, 
    { "COMPOUND_TEXT",   0, TARGET_COMPOUND_TEXT }
};

GtkTargetEntry rtf_targets[] = {
    { "application/rtf",   0, TARGET_RTF },
    { "application/x-rtf", 0, TARGET_RTF },
    { "text/rtf",          0, TARGET_RTF },
    { "text/richtext",     0, TARGET_RTF },
    { "STRING",            0, TARGET_STRING },
    { "TEXT",              0, TARGET_TEXT }
};

static int n_text = sizeof text_targets / sizeof text_targets[0];
static int n_rtf = sizeof rtf_targets / sizeof rtf_targets[0];

static void gretl_clipboard_set (int copycode);

static void gretl_clipboard_free (void)
{
    free(clipboard_buf);
    clipboard_buf = NULL;
}

int buf_to_clipboard (const char *buf)
{
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();
    clipboard_buf = gretl_strdup(buf);
    if (clipboard_buf == NULL) {
	err = 1;
    } else {
	gretl_clipboard_set(GRETL_FORMAT_TXT);
    }

    return err;
}

#if CLIPDEBUG

static const char *fmt_label (int f)
{
    if (f == TARGET_UTF8_STRING) {
	return "TARGET_UTF8_STRING";
    } else if (f == TARGET_STRING) {
	return "TARGET_STRING";
    } else if (f == TARGET_TEXT) {
	return "TARGET_STRING";
    } else if (f == TARGET_COMPOUND_TEXT) {
	return "TARGET_COMPOUND_STRING";
    } else if (f == TARGET_RTF) {
	return "TARGET_RTF";
    } else {
	return "unknown";
    }
}

#endif  

static void gretl_clipboard_get (GtkClipboard *clip,
				 GtkSelectionData *selection_data,
				 guint info,
				 gpointer p)
{
    gchar *str = clipboard_buf; /* global */

#if CLIPDEBUG
    fprintf(stderr, "gretl_clipboard_get: info = %d (%s)\n", 
	    (int) info, fmt_label(info));
#endif   

    if (str == NULL || *str == '\0') {
	return;
    }

    if (info != TARGET_UTF8_STRING) {
	/* remove any Unicode minuses */
	strip_unicode_minus(str);
    }

    if (info == TARGET_RTF) {
	gtk_selection_data_set(selection_data,
			       GDK_SELECTION_TYPE_STRING,
			       8, (guchar *) str, 
			       strlen(str));
    } else {
	gtk_selection_data_set_text(selection_data, str, -1);
    }
}

static void gretl_clipboard_clear (GtkClipboard *clip, gpointer p)
{
    gretl_clipboard_free();
}

#ifdef OS_OSX

static void pasteboard_set (int fmt)
{
    FILE *fp = popen("/usr/bin/pbcopy", "w");

    if (fp == NULL) {
	errbox("Couldn't open pbcopy");
    } else {
	fputs(clipboard_buf, fp);
	pclose(fp);
    }
}

#endif

static void gretl_clipboard_set (int fmt)
{
    static GtkClipboard *clip;
    GtkTargetEntry *targs;
    gint n_targs;

#ifdef OS_OSX
    if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) {
	pasteboard_set(fmt);
	return;
    }
#endif

    if (clip == NULL) {
	clip = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
    }

    if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) {
	targs = rtf_targets;
	n_targs = n_rtf;
    } else {
	targs = text_targets;
	n_targs = n_text;
    }

#if CLIPDEBUG
    fprintf(stderr, "gretl_clipboard_set: fmt = %d\n", fmt);
#endif

    if (!gtk_clipboard_set_with_owner(clip, targs, n_targs,
				      gretl_clipboard_get,
				      gretl_clipboard_clear,
				      G_OBJECT(mdata->main))) {
	fprintf(stderr, "Failed to initialize clipboard\n");
    }
}

/* We call this when a buffer to be copied as RTF contains
   non-ASCII characters but validates as UTF-8. Note that the 
   UTF-8 minus sign does not recode correctly into Windows 
   codepages, so if it's present we have to remove it first, 
   then see if there's anything non-ASCII left to handle.
*/

static char *strip_utf8 (char *buf)
{
    strip_unicode_minus(buf);

    if (string_is_utf8((const unsigned char *) buf)) {
#if CLIPDEBUG
	fprintf(stderr, "strip_utf8: stage 2, recoding\n");
#endif
	return utf8_to_cp(buf);
    } else {
	return buf;
    }
}

/* note: there's a Windows-specific counterpart to this
   in gretlwin32.c
*/

int prn_to_clipboard (PRN *prn, int fmt)
{
    char *buf = gretl_print_steal_buffer(prn);
    const char *cset = NULL;
    gchar *trbuf = NULL;
    int utf8_ok = 1;
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();

    if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) {
	/* UTF-8 will never work in RTF */
	utf8_ok = 0;
    } else if (fmt == GRETL_FORMAT_TXT && !g_get_charset(&cset)) {
	utf8_ok = 0;
    }

#if CLIPDEBUG
    fprintf(stderr, "prn_to_clipboard, fmt = %d, utf8_ok = %d\n", fmt, utf8_ok);
#endif

    if (!utf8_ok && string_is_utf8((const unsigned char *) buf)) {
	trbuf = strip_utf8(buf);
    } else {
	trbuf = buf;
    }

    if (fmt == GRETL_FORMAT_RTF_TXT) {
	clipboard_buf = dosify_buffer(trbuf, fmt);
    } else {
	clipboard_buf = gretl_strdup(trbuf);
    }

    if (trbuf != buf) {
	g_free(trbuf);
    }

    if (clipboard_buf == NULL) {
	err = 1;
    } else {
	gretl_clipboard_set(fmt);
    }

    return err;
}


