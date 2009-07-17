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

static gchar *clipboard_buf; 

GtkTargetEntry basic_targets[] = {
    { "UTF8_STRING",     0, TARGET_UTF8_STRING },
    { "STRING",          0, TARGET_STRING },
    { "TEXT",            0, TARGET_TEXT }, 
    { "COMPOUND_TEXT",   0, TARGET_COMPOUND_TEXT }
};

GtkTargetEntry rtf_targets[] = {
    { "application/rtf", 0, TARGET_RTF },   
    { "text/rtf",        0, TARGET_RTF }
};

static int n_basic = sizeof basic_targets / sizeof basic_targets[0];
static int n_rtf = sizeof rtf_targets / sizeof rtf_targets[0];

static void gretl_clipboard_set (int copycode);

#define CLIPDEBUG 0

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

#define UTF_MINUS(u,i) (u[i] == 0xE2 && u[i+1] == 0x88 && u[i+2] == 0x92)

static gchar *minus_check (gchar *s)
{
    guchar *u = (guchar *) s;
    gchar *ret = s;
    int i, n = strlen(s);
    int got_minus = 0;

    for (i=0; i<n-3; i++) {
	if (UTF_MINUS(u, i)) {
	    got_minus = 1;
	    break;
	}
    }

    if (got_minus) {
	/* convert U+2212 into plain ASCII dash */
	int j = 0;

	ret = calloc(n, 1);
	if (ret != NULL) {
	    for (i=0; i<n; i++) {
		if (i < n - 3 && UTF_MINUS(u, i)) {
		    ret[j++] = '-';
		    i += 2;
		} else {
		    ret[j++] = u[i];
		}
	    }
	} else {
	    ret = s;
	}
    } 

    return ret;
}

static void gretl_clipboard_get (GtkClipboard *clip,
				 GtkSelectionData *selection_data,
				 guint info,
				 gpointer p)
{
    gchar *str = clipboard_buf; /* global */

#if CLIPDEBUG
    fprintf(stderr, "info = %d\n", (int) info);
    if (info == TARGET_UTF8_STRING) {
	fprintf(stderr, " = TARGET_UTF8_STRING\n");
    } else if (info == TARGET_STRING) {
	fprintf(stderr, " = TARGET_STRING\n");
    } else if (info == TARGET_TEXT) {
	fprintf(stderr, " = TARGET_STRING\n");
    } else if (info == TARGET_COMPOUND_TEXT) {
	fprintf(stderr, " = TARGET_COMPOUND_STRING\n");
    } else if (info == TARGET_RTF) {
	fprintf(stderr, " = TARGET_RTF\n");
    }
#endif   

    if (str == NULL || *str == '\0') {
	return;
    }

    if (info != TARGET_UTF8_STRING) {
	/* need to remove any UTF-8 minuses? */
	str = minus_check(str);
    }

    if (info == TARGET_RTF) {
	gtk_selection_data_set(selection_data,
			       GDK_SELECTION_TYPE_STRING,
			       8 * sizeof(gchar), 
			       (guchar *) str, 
			       strlen(str));
    } else {
	gtk_selection_data_set_text(selection_data, str, -1);
    }

    if (str != clipboard_buf) {
	g_free(str);
    }
}

static void gretl_clipboard_clear (GtkClipboard *clip, gpointer p)
{
    gretl_clipboard_free();
}

static void gretl_clipboard_set (int fmt)
{
    static GtkClipboard *clip;
    GtkTargetEntry *targs;
    gint n_targs;

    if (clip == NULL) {
	clip = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
    }

    if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) {
	targs = rtf_targets;
	n_targs = n_rtf;
    } else {
	targs = basic_targets;
	n_targs = n_basic;
    }

    if (!gtk_clipboard_set_with_owner(clip, targs, n_targs,
				      gretl_clipboard_get,
				      gretl_clipboard_clear,
				      G_OBJECT(mdata->main))) {
	fprintf(stderr, "Failed to initialize clipboard\n");
    }
}

int prn_to_clipboard (PRN *prn, int fmt)
{
    const char *buf = gretl_print_get_buffer(prn);
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();

    if (fmt == GRETL_FORMAT_TXT) {
	/* recode (only) if needed */
	clipboard_buf = my_locale_from_utf8(buf);
    } else if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) { 
	/* RTF: ensure that we're not in UTF-8 */
	if (string_is_utf8((const unsigned char *) buf)) {
	    gchar *trbuf = utf8_to_cp(buf);

	    if (trbuf != NULL) {
		if (fmt == GRETL_FORMAT_RTF_TXT) {
		    clipboard_buf = dosify_buffer(trbuf, fmt);
		} else {
		    clipboard_buf = gretl_strdup(trbuf);
		}
		g_free(trbuf);
	    }
	} else if (fmt == GRETL_FORMAT_RTF_TXT) {
	    clipboard_buf = dosify_buffer(buf, fmt);
	} else {
	    clipboard_buf = gretl_strdup(buf);
	}
    } else {
	/* TeX, CSV */
	clipboard_buf = gretl_strdup(buf);
    }

    if (clipboard_buf == NULL) {
	err = 1;
    } else {
	gretl_clipboard_set(fmt);
    }

    return err;
}


