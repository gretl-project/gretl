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

gchar *clipboard_buf; 

GtkTargetEntry basic_targets[] = {
    { "STRING",          0, TARGET_STRING },
    { "TEXT",            0, TARGET_TEXT }, 
    { "COMPOUND_TEXT",   0, TARGET_COMPOUND_TEXT }
};

GtkTargetEntry extended_targets[] = {
    { "STRING",          0, TARGET_STRING },
    { "TEXT",            0, TARGET_TEXT }, 
    { "COMPOUND_TEXT",   0, TARGET_COMPOUND_TEXT },
    { "text/rtf",        0, TARGET_RTF },
    { "application/rtf", 0, TARGET_RTF }
};

enum {
    n_basic_targets    = 3,
    n_extended_targets = 5
};

#define CLIPDEBUG

void gretl_clipboard_free (void)
{
    free(clipboard_buf);
    clipboard_buf = NULL;
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

static void gretl_clip_init (int copycode)
{
    static GtkClipboard *clip;
    GtkTargetEntry *targets;
    gint n_targets;

    if (copycode == COPY_RTF || copycode == COPY_TEXT_AS_RTF) {
	targets = extended_targets;
	n_targets = n_extended_targets;
    } else {
	targets = basic_targets;
	n_targets = n_basic_targets;
    }

    if (clip == NULL) {
	clip = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
    }

    if (!gtk_clipboard_set_with_owner(clip,
				      targets, n_targets,
				      gretl_clipboard_get,
				      gretl_clipboard_clear,
				      G_OBJECT(mdata->w))) {
	fprintf(stderr, "Failed to initialize clipboard\n");
    }
}

void gretl_clipboard_set (int copycode)
{
    gretl_clip_init(copycode);

    gtk_selection_owner_set(mdata->w,
			    GDK_SELECTION_PRIMARY, 
			    GDK_CURRENT_TIME);
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

static void gretl_clip_init (int copycode)
{
    GdkAtom clipboard_atom = GDK_NONE;
    GtkTargetEntry *targets;
    gint n_targets;

    if (1 || copycode == COPY_RTF || copycode == COPY_TEXT_AS_RTF) {
	targets = extended_targets;
	n_targets = n_extended_targets;
    } else {
	targets = basic_targets;
	n_targets = n_basic_targets;
    }

    clipboard_atom = gdk_atom_intern("CLIPBOARD", FALSE);

    gtk_selection_add_targets(mdata->w, GDK_SELECTION_PRIMARY,
			      targets, n_targets);
    gtk_selection_add_targets(mdata->w, clipboard_atom,
			      targets, n_targets);
    gtk_signal_connect(GTK_OBJECT(mdata->w), "selection_get",
		       GTK_SIGNAL_FUNC(gretl_clipboard_get), NULL);    
}

void gretl_clipboard_set (int copycode)
{
    static int initted;

    if (!initted) {
	gretl_clip_init(copycode);
	initted = 1;
    }

    gtk_selection_owner_set(mdata->w,
			    GDK_SELECTION_PRIMARY, 
			    GDK_CURRENT_TIME);
}

#endif
