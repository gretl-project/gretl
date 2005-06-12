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
#include "textutil.h"

static gchar *clipboard_buf; 

GtkTargetEntry basic_targets[] = {
    { "STRING",          0, TARGET_STRING },
    { "TEXT",            0, TARGET_TEXT }, 
    { "COMPOUND_TEXT",   0, TARGET_COMPOUND_TEXT }
};

GtkTargetEntry full_targets[] = {
    { "STRING",          0, TARGET_STRING },
    { "TEXT",            0, TARGET_TEXT }, 
    { "COMPOUND_TEXT",   0, TARGET_COMPOUND_TEXT },
    { "text/rtf",        0, TARGET_RTF },
    { "application/rtf", 0, TARGET_RTF }
};

#ifndef OLD_GTK
static int n_basic = sizeof basic_targets / sizeof basic_targets[0];
#endif
static int n_full = sizeof full_targets / sizeof full_targets[0];

static void gretl_clipboard_set (int copycode);

#undef CLIPDEBUG

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

#ifndef OLD_GTK

static void gretl_clipboard_get (GtkClipboard *clip,
				 GtkSelectionData *selection_data,
				 guint info,
				 gpointer p)
{
    gchar *str;
    gint length;

#ifdef CLIPDEBUG
    fprintf(stderr, "info = %d\n", (int) info);
    if (info == TARGET_STRING) {
	fprintf(stderr, " = TARGET_STRING\n");
    } else if (info == TARGET_TEXT) {
	fprintf(stderr, " = TARGET_STRING\n");
    } else if (info == TARGET_COMPOUND_TEXT) {
	fprintf(stderr, " = TARGET_COMPOUND_STRING\n");
    } else if (info == TARGET_RTF) {
	fprintf(stderr, " = TARGET_RTF\n");
    }
#endif    

    str = clipboard_buf; /* global */
    if (str == NULL) {
	return;
    }

    length = strlen(str);

    if (info == TARGET_STRING || info == TARGET_RTF) {
	gtk_selection_data_set (selection_data,
				GDK_SELECTION_TYPE_STRING,
				8 * sizeof(gchar), 
				(guchar *) str, 
				length);
    } else if (info == TARGET_TEXT || info == TARGET_COMPOUND_TEXT) {
	guchar *text;
	gchar c;
	GdkAtom seltype;
	gint format;
	gint new_length;

	c = str[length];
	str[length] = '\0';
	gdk_string_to_compound_text(str, &seltype, &format, 
				    &text, &new_length);
	gtk_selection_data_set(selection_data, seltype, format, 
			       text, new_length);
	gdk_free_compound_text(text);
	str[length] = c;
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
	targs = full_targets;
	n_targs = n_full;
    } else {
	targs = basic_targets;
	n_targs = n_basic;
    }

    if (!gtk_clipboard_set_with_owner(clip, targs, n_targs,
				      gretl_clipboard_get,
				      gretl_clipboard_clear,
				      G_OBJECT(mdata->w))) {
	fprintf(stderr, "Failed to initialize clipboard\n");
    }
}

#else /* gtk-1.2 version */

static gint 
gretl_clipboard_get (GtkWidget *widget,
		     GtkSelectionData *selection_data,
		     guint info,
		     guint time)
{
    gchar *str;
    gint length;

    str = clipboard_buf;
    if (str == NULL) {
	return TRUE;
    }

    length = strlen(str);

#ifdef CLIPDEBUG
    fprintf(stderr, "info = %d\n", (int) info);
    if (info == TARGET_STRING) {
	fprintf(stderr, " = TARGET_STRING\n");
    } else if (info == TARGET_TEXT) {
	fprintf(stderr, " = TARGET_STRING\n");
    } else if (info == TARGET_COMPOUND_TEXT) {
	fprintf(stderr, " = TARGET_COMPOUND_STRING\n");
    } else if (info == TARGET_RTF) {
	fprintf(stderr, " = TARGET_RTF\n");
    }
#endif    
  
    if (info == TARGET_STRING || info == TARGET_RTF) {
	gtk_selection_data_set (selection_data,
				GDK_SELECTION_TYPE_STRING,
				8 * sizeof(gchar), 
				(guchar *) str, 
				length);
    } else if (info == TARGET_TEXT || info == TARGET_COMPOUND_TEXT) {
	guchar *text;
	gchar c;
	GdkAtom encoding;
	gint format;
	gint new_length;

	c = str[length];
	str[length] = '\0';
	gdk_string_to_compound_text(str, &encoding, &format, 
				    &text, &new_length);
	gtk_selection_data_set(selection_data, encoding, format, 
			       text, new_length);
	gdk_free_compound_text(text);
	str[length] = c;
    }

    return TRUE;
}

static void gretl_clip_init (int copycode, GdkAtom clipatom)
{
    gtk_selection_add_targets(mdata->w, clipatom, full_targets, n_full);
    gtk_signal_connect(GTK_OBJECT(mdata->w), "selection_get",
		       GTK_SIGNAL_FUNC(gretl_clipboard_get), NULL); 
}

static void gretl_clipboard_set (int copycode)
{
    GdkAtom clipatom = GDK_NONE;

    if (clipatom == GDK_NONE) {
	clipatom = gdk_atom_intern("CLIPBOARD", FALSE);
	gretl_clip_init(copycode, clipatom);
    }

    gtk_selection_owner_set(mdata->w, clipatom, GDK_CURRENT_TIME);      
}

#endif

#if defined(ENABLE_NLS) && !defined(OLD_GTK)

int prn_to_clipboard (PRN *prn, int fmt)
{
    const char *buf = gretl_print_get_buffer(prn);
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();

    if (fmt == GRETL_FORMAT_TXT || fmt == GRETL_FORMAT_RTF_TXT) { 
	/* need to convert from utf8 */
	gchar *trbuf = my_locale_from_utf8(buf);

	if (trbuf == NULL) {
	    err = 1;
	} else if (fmt == GRETL_FORMAT_TXT) {
	    clipboard_buf = gretl_strdup(trbuf);
	    if (clipboard_buf == NULL) {
		err = 1;
	    }
	    g_free(trbuf);
	} else if (fmt == GRETL_FORMAT_RTF_TXT) {
	    clipboard_buf = dosify_buffer(trbuf, fmt);
	    if (clipboard_buf == NULL) {
		err = 1;
	    }
	    g_free(trbuf);
	}
    } else { 
	/* copying TeX, RTF or CSV */
	clipboard_buf = gretl_strdup(buf);
	err = (clipboard_buf == NULL);
    }

    if (!err) {
	gretl_clipboard_set(fmt);
    }

    return err;
}

#else /* plain GTK, no utf-8 */

int prn_to_clipboard (PRN *prn, int copycode)
{
    const char *buf;
    int err = 0;

    buf = gretl_print_get_buffer(prn);

    if (buf == NULL || *buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();
    clipboard_buf = gretl_strdup(buf);
    err = (clipboard_buf == NULL);

    if (!err) {
	gretl_clipboard_set(copycode);
    }

    return err;
}

#endif /* switch for prn_to_clipboard */


