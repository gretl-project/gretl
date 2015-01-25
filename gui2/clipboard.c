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

#define CLIPDEBUG 2

static gchar *clipboard_buf;
static gsize clipboard_bytes;

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

GtkTargetEntry svg_targets[] = {
    { "application/svg+xml", 0, TARGET_SVG },
    { "image/svg+xml",       0, TARGET_SVG },
    { "image/svg",           0, TARGET_SVG }
};

GtkTargetEntry emf_targets[] = {
    { "application/emf",   0, TARGET_EMF },
    { "application/x-emf", 0, TARGET_EMF },
    { "image/x-emf",       0, TARGET_EMF }
};

static int n_text = sizeof text_targets / sizeof text_targets[0];
static int n_rtf = sizeof rtf_targets / sizeof rtf_targets[0];
static int n_svg = sizeof svg_targets / sizeof svg_targets[0];
static int n_emf = sizeof emf_targets / sizeof emf_targets[0];

static void gretl_clipboard_set (int copycode, int svg, int emf);

static void gretl_clipboard_free (void)
{
    free(clipboard_buf);
    clipboard_buf = NULL;
    clipboard_bytes = 0;
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
	gretl_clipboard_set(GRETL_FORMAT_TXT, 0, 0);
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
    } else if (f == TARGET_SVG) {
	return "TARGET_SVG";
    } else if (f == TARGET_EMF) {
	return "TARGET_EMF";
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
#if CLIPDEBUG
    fprintf(stderr, "gretl_clipboard_get: info = %d (%s)\n", 
	    (int) info, fmt_label(info));
#endif   

    if (clipboard_buf == NULL || *clipboard_buf == '\0') {
	return;
    }

    if (info != TARGET_UTF8_STRING &&
	info != TARGET_SVG && info != TARGET_EMF) {
	/* remove any Unicode minuses (??) */
	strip_unicode_minus(clipboard_buf);
    }

    if (info == TARGET_RTF || info == TARGET_SVG) {
	gtk_selection_data_set(selection_data,
			       GDK_SELECTION_TYPE_STRING,
			       8, (guchar *) clipboard_buf, 
			       strlen(clipboard_buf));
    } else if (info == TARGET_EMF) {
	gtk_selection_data_set(selection_data,
			       GDK_SELECTION_TYPE_STRING,
			       8, (guchar *) clipboard_buf, 
			       clipboard_bytes);
    } else {
	gtk_selection_data_set_text(selection_data, clipboard_buf, -1);
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

static void gretl_clipboard_set (int fmt, int svg, int emf)
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

    if (svg) {
	targs = svg_targets;
	n_targs = n_svg;
    } else if (emf) {
	targs = emf_targets;
	n_targs = n_emf;
    } else if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) {
	targs = rtf_targets;
	n_targs = n_rtf;
    } else {
	targs = text_targets;
	n_targs = n_text;
    }

#if CLIPDEBUG
    fprintf(stderr, "gretl_clipboard_set: fmt = %d, svg = %d\n", fmt, svg);
#endif

    if (!gtk_clipboard_set_with_owner(clip, targs, n_targs,
				      gretl_clipboard_get,
				      gretl_clipboard_clear,
				      G_OBJECT(mdata->main))) {
	fprintf(stderr, "Failed to initialize clipboard\n");
    }
}

/* note: there's a Windows-specific counterpart to this
   in gretlwin32.c
*/

int prn_to_clipboard (PRN *prn, int fmt)
{
    char *buf = gretl_print_steal_buffer(prn);
    char *modbuf = NULL;
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();

    err = maybe_post_process_buffer(buf, fmt, W_COPY, &modbuf);

    if (!err) {
	if (modbuf != NULL) {
	    clipboard_buf = modbuf;
	} else {
	    clipboard_buf = buf;
	}
	gretl_clipboard_set(fmt, 0, 0);
    }

    if (buf != clipboard_buf) {
	free(buf);
    }

    return err;
}

int svg_to_clipboard (const char *fname)
{
    gchar *buf = NULL;
    gsize sz = 0;
    
    gretl_file_get_contents(fname, &buf, &sz);

    fprintf(stderr, "clip: got buf at %p, size %d\n",
	    (void *) buf, (int) sz);

    if (buf != NULL && *buf != '\0') {
	gretl_clipboard_free();
	clipboard_buf = buf;
	clipboard_bytes = sz;
	gretl_clipboard_set(0, 1, 0);
    }

    return 0;
}


