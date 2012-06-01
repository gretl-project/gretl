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

/*
 * Modified by the GTK+ Team and others 1997-2000.  See the AUTHORS
 * file for a list of people on the GTK+ Team.  See the ChangeLog
 * files for a list of changes.  These files are distributed with
 * GTK+ at ftp://ftp.gtk.org/pub/gtk/. 
 */

/*
 * Modified to allow filtering of fonts by Allin Cottrell, 2002
 */

#ifndef __GTK_FONTSELHACK_H__
#define __GTK_FONTSELHACK_H__

#include <gdk/gdk.h>
#include <glib-object.h>

#define GTK_TYPE_FONTSEL_HACK_DIALOG  (gtk_fontsel_hack_dialog_get_type ())
#define GTK_FONTSEL_HACK_DIALOG(obj)  (G_TYPE_CHECK_INSTANCE_CAST((obj), GTK_TYPE_FONTSEL_HACK_DIALOG, \
				       GtkFontselHackDialog))

typedef struct _GtkFontselHack GtkFontselHack;
typedef struct _GtkFontselHackClass GtkFontselHackClass;

typedef struct _GtkFontselHackDialog GtkFontselHackDialog;
typedef struct _GtkFontselHackDialogClass GtkFontselHackDialogClass;

typedef enum {
    FONT_HACK_NONE,
    FONT_HACK_LATIN,
    FONT_HACK_LATIN_MONO
} FontFilterType;

GType	   gtk_fontsel_hack_get_type		(void) G_GNUC_CONST;
GtkWidget* gtk_fontsel_hack_new		        (void);
gchar*	   gtk_fontsel_hack_get_font_name	(GtkFontselHack *fontsel);

gboolean              gtk_fontsel_hack_set_font_name    (GtkFontselHack *fontsel,
							 const gchar      *fontname);
G_CONST_RETURN gchar* gtk_fontsel_hack_get_preview_text (GtkFontselHack *fontsel);
void                  gtk_fontsel_hack_set_preview_text (GtkFontselHack *fontsel,
							 const gchar      *text);

gint gtk_fontsel_hack_get_filter (GtkFontselHack *fontsel);
void gtk_fontsel_hack_set_filter (GtkFontselHack *fontsel,
				  FontFilterType filter);

GType gtk_fontsel_hack_dialog_get_type (void) G_GNUC_CONST;

GtkWidget* gtk_fontsel_hack_dialog_new (const gchar *title);

gchar* gtk_fontsel_hack_dialog_get_font_name(GtkFontselHackDialog *fsd);

gboolean gtk_fontsel_hack_dialog_set_font_name(GtkFontselHackDialog *fsd,
					       const gchar *fontname);

G_CONST_RETURN gchar* gtk_fontsel_hack_dialog_get_preview_text (GtkFontselHackDialog *fsd);

/* This sets the text in the preview entry. It will be copied by the entry. */
void gtk_fontsel_hack_dialog_set_preview_text(GtkFontselHackDialog *fsd,
					      const gchar *text);

/* This returns the font filter currently in force. */
gint gtk_fontsel_hack_dialog_get_filter(GtkFontselHackDialog *fsd);

/* This sets the font filter for the dialog. */
void gtk_fontsel_hack_dialog_set_filter (GtkFontselHackDialog *fsd,
					 FontFilterType filter);

GtkWidget *gtk_fontsel_hack_dialog_ok_button(GtkWidget *fsd);

GtkWidget *gtk_fontsel_hack_dialog_cancel_button(GtkWidget *fsd);


#endif /* __GTK_FONTSELHACK_H__ */
