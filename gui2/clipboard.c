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

static int n_basic = sizeof basic_targets / sizeof basic_targets[0];
static int n_full = sizeof full_targets / sizeof full_targets[0];

static void gretl_clipboard_set (int copycode);

#define CLIPDEBUG

static void gretl_clipboard_free (void)
{
    free(clipboard_buf);
    clipboard_buf = NULL;
}

static int copy_to_clipboard_buf (const char *buf)
{
    size_t len = strlen(buf);

    clipboard_buf = mymalloc(len + 1);
    if (clipboard_buf == NULL) {
	return 1;
    }

    memcpy(clipboard_buf, buf, len + 1);

    return 0;
}

int buf_to_clipboard (const char *buf)
{
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();

    err = copy_to_clipboard_buf(buf);

    if (!err) {
	gretl_clipboard_set(COPY_TEXT);
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

static void gretl_clipboard_set (int copycode)
{
    static GtkClipboard *clip;
    GtkTargetEntry *targs;
    gint n_targs;

    if (clip == NULL) {
	clip = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
    }

    if (copycode == COPY_RTF || copycode == COPY_TEXT_AS_RTF) {
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
    GtkTargetEntry *targs;
    gint n_targs;

    if (1 || copycode == COPY_RTF || copycode == COPY_TEXT_AS_RTF) {
	targs = full_targets;
	n_targs = n_full;
    } else {
	targs = basic_targets;
	n_targs = n_basic;
    }

    gtk_selection_add_targets(mdata->w, clipatom, targs, n_targs);
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

int prn_to_clipboard (PRN *prn, int copycode)
{
    if (prn->buf == NULL || *prn->buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();

    if (copycode == COPY_TEXT || copycode == COPY_TEXT_AS_RTF) { 
	/* need to convert from utf8 */
	gchar *trbuf = my_locale_from_utf8(prn->buf);
	
	if (trbuf != NULL) {
	    size_t len = strlen(trbuf);

	    if (copycode == COPY_TEXT_AS_RTF) {
		clipboard_buf = dosify_buffer(trbuf, copycode);
	    } else {
		clipboard_buf = mymalloc(len + 1);
	    }
	    if (clipboard_buf == NULL) {
		g_free(trbuf);
		return 1;
	    }
	    if (copycode != COPY_TEXT_AS_RTF) {
		memcpy(clipboard_buf, trbuf, len + 1);
	    }
	    g_free(trbuf);
	}
    } else { /* copying TeX, RTF or CSV */
	if (copy_to_clipboard_buf(prn->buf)) {
	    return 1;
	}
    }

    gretl_clipboard_set(copycode);

    return 0;
}

#else /* plain GTK, no NLS */

int prn_to_clipboard (PRN *prn, int copycode)
{
    int err = 0;

    if (prn->buf == NULL || *prn->buf == '\0') {
	return 0;
    }

    gretl_clipboard_free();

    err = copy_to_clipboard_buf(prn->buf);

    if (!err) {
	gretl_clipboard_set(copycode);
    }

    return err;
}

#endif /* switch for prn_to_clipboard */


