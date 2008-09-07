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
#include <gtk/gtkdialog.h>
#include <gtk/gtkvbox.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define GTK_TYPE_FONT_SELECTION_HACK              (gtk_font_selection_hack_get_type ())
#define GTK_FONT_SELECTION_HACK(obj)              (GTK_CHECK_CAST ((obj), GTK_TYPE_FONT_SELECTION_HACK, GtkFontSelectionHack))
#define GTK_FONT_SELECTION_HACK_CLASS(klass)      (GTK_CHECK_CLASS_CAST ((klass), GTK_TYPE_FONT_SELECTION_HACK, GtkFontSelectionHackClass))
#define GTK_IS_FONT_SELECTION_HACK(obj)           (GTK_CHECK_TYPE ((obj), GTK_TYPE_FONT_SELECTION_HACK))
#define GTK_IS_FONT_SELECTION_HACK_CLASS(klass)   (GTK_CHECK_CLASS_TYPE ((klass), GTK_TYPE_FONT_SELECTION_HACK))
#define GTK_FONT_SELECTION_HACK_GET_CLASS(obj)    (GTK_CHECK_GET_CLASS ((obj), GTK_TYPE_FONT_SELECTION_HACK, GtkFontSelectionHackClass))


#define GTK_TYPE_FONT_SELECTION_HACK_DIALOG              (gtk_font_selection_hack_dialog_get_type ())
#define GTK_FONT_SELECTION_HACK_DIALOG(obj)              (GTK_CHECK_CAST ((obj), GTK_TYPE_FONT_SELECTION_HACK_DIALOG, GtkFontSelectionHackDialog))
#define GTK_FONT_SELECTION_HACK_DIALOG_CLASS(klass)      (GTK_CHECK_CLASS_CAST ((klass), GTK_TYPE_FONT_SELECTION_HACK_DIALOG, GtkFontSelectionHackDialogClass))
#define GTK_IS_FONT_SELECTION_HACK_DIALOG(obj)           (GTK_CHECK_TYPE ((obj), GTK_TYPE_FONT_SELECTION_HACK_DIALOG))
#define GTK_IS_FONT_SELECTION_HACK_DIALOG_CLASS(klass)   (GTK_CHECK_CLASS_TYPE ((klass), GTK_TYPE_FONT_SELECTION_HACK_DIALOG))
#define GTK_FONT_SELECTION_HACK_DIALOG_GET_CLASS(obj)    (GTK_CHECK_GET_CLASS ((obj), GTK_TYPE_FONT_SELECTION_HACK_DIALOG, GtkFontSelectionHackDialogClass))

#define GTK_TYPE_FONT_FILTER (gtk_font_filter_get_type())


  typedef struct _GtkFontSelectionHack	     GtkFontSelectionHack;
  typedef struct _GtkFontSelectionHackClass	     GtkFontSelectionHackClass;

  typedef struct _GtkFontSelectionHackDialog	     GtkFontSelectionHackDialog;
  typedef struct _GtkFontSelectionHackDialogClass  GtkFontSelectionHackDialogClass;

  typedef enum {
    GTK_FONT_HACK_NONE,
    GTK_FONT_HACK_LATIN,
    GTK_FONT_HACK_LATIN_MONO
  } GtkFontFilterType;

  struct _GtkFontSelectionHack
  {
    GtkVBox parent_instance;
  
    GtkWidget *font_entry;
    GtkWidget *family_list;
    GtkWidget *font_style_entry;
    GtkWidget *face_list;
    GtkWidget *size_entry;
    GtkWidget *size_list;
    GtkWidget *pixels_button;
    GtkWidget *points_button;
    GtkWidget *filter_button;
    GtkWidget *preview_entry;

    PangoFontFamily *family;	/* Current family */
    PangoFontFace *face;		/* Current face */
  
    gint size;
    GtkFontFilterType filter;
  
    GdkFont *font;		/* Cache for gdk_font_selection_hack_get_font, so we can preserve
				 * refcounting behavior
				 */
  };

  struct _GtkFontSelectionHackClass
  {
    GtkVBoxClass parent_class;

    /* Padding for future expansion */
    void (*_gtk_reserved1) (void);
    void (*_gtk_reserved2) (void);
    void (*_gtk_reserved3) (void);
    void (*_gtk_reserved4) (void);
  };


  struct _GtkFontSelectionHackDialog
  {
    GtkDialog parent_instance;
  
    GtkWidget *fontsel;
  
    GtkWidget *main_vbox;
    GtkWidget *action_area;
    GtkWidget *ok_button;
    /* The 'Apply' button is not shown by default but you can show/hide it. */
    GtkWidget *apply_button;
    GtkWidget *cancel_button;
  
    /* If the user changes the width of the dialog, we turn auto-shrink off. */
    gint dialog_width;
    gboolean auto_resize;
  };

  struct _GtkFontSelectionHackDialogClass
  {
    GtkDialogClass parent_class;

    /* Padding for future expansion */
    void (*_gtk_reserved1) (void);
    void (*_gtk_reserved2) (void);
    void (*_gtk_reserved3) (void);
    void (*_gtk_reserved4) (void);
  };



  /*****************************************************************************
   * GtkFontSelectionHack functions.
   *   see the comments in the GtkFontSelectionHackDialog functions.
   *****************************************************************************/

  GtkType	   gtk_font_selection_hack_get_type		(void) G_GNUC_CONST;
  GtkWidget* gtk_font_selection_hack_new		(void);
  gchar*	   gtk_font_selection_hack_get_font_name	(GtkFontSelectionHack *fontsel);

  gboolean              gtk_font_selection_hack_set_font_name    (GtkFontSelectionHack *fontsel,
								  const gchar      *fontname);
  G_CONST_RETURN gchar* gtk_font_selection_hack_get_preview_text (GtkFontSelectionHack *fontsel);
  void                  gtk_font_selection_hack_set_preview_text (GtkFontSelectionHack *fontsel,
								  const gchar      *text);
  gint                  gtk_font_selection_hack_get_filter        (GtkFontSelectionHack *fontsel);
  void                  gtk_font_selection_hack_set_filter        (GtkFontSelectionHack *fontsel,
								   GtkFontFilterType      filter);


  /*****************************************************************************
   * GtkFontSelectionHackDialog functions.
   *   most of these functions simply call the corresponding function in the
   *   GtkFontSelectionHack.
   *****************************************************************************/

  GtkType	   gtk_font_selection_hack_dialog_get_type	(void) G_GNUC_CONST;
    GtkWidget* gtk_font_selection_hack_dialog_new	(const gchar	  *title);

  /* This returns the X Logical Font Description fontname, or NULL if no font
     is selected. Note that there is a slight possibility that the font might not
     have been loaded OK. You should call gtk_font_selection_hack_dialog_get_font()
     to see if it has been loaded OK.
     You should g_free() the returned font name after you're done with it. */
  gchar*	 gtk_font_selection_hack_dialog_get_font_name    (GtkFontSelectionHackDialog *fsd);

  /* This sets the currently displayed font. It should be a valid X Logical
     Font Description font name (anything else will be ignored), e.g.
     "-adobe-courier-bold-o-normal--25-*-*-*-*-*-*-*" 
     It returns TRUE on success. */
  gboolean gtk_font_selection_hack_dialog_set_font_name    (GtkFontSelectionHackDialog *fsd,
							    const gchar	*fontname);

  /* This returns the text in the preview entry. You should copy the returned
     text if you need it. */
  G_CONST_RETURN gchar* gtk_font_selection_hack_dialog_get_preview_text (GtkFontSelectionHackDialog *fsd);

  /* This sets the text in the preview entry. It will be copied by the entry,
     so there's no need to g_strdup() it first. */
  void	 gtk_font_selection_hack_dialog_set_preview_text (GtkFontSelectionHackDialog *fsd,
							  const gchar	    *text);
  /* This returns the font filter currently in force. */
  gint   gtk_font_selection_hack_dialog_get_filter (GtkFontSelectionHackDialog *fsd);

  /* This sets the font filter for the dialog. */
  void   gtk_font_selection_hack_dialog_set_filter (GtkFontSelectionHackDialog *fsd,
						    GtkFontFilterType      filter);

  int font_has_minus (PangoFontDescription *desc);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __GTK_FONTSELHACK_H__ */
