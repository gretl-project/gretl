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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef DLGUTILS_H
#define DLGUTILS_H

GtkWidget *gretl_dialog_new (const char *title);

void set_dialog_border_widths (GtkWidget *dlg);

void dialog_set_no_resize (GtkWidget *w);

GtkWidget *context_help_button (GtkWidget *hbox, int cmdcode);

GtkWidget *cancel_delete_button (GtkWidget *hbox, GtkWidget *targ);

GtkWidget *cancel_options_button (GtkWidget *hbox, GtkWidget *targ,
				  int *opt);

GtkWidget *ok_button (GtkWidget *hbox);

GtkWidget *next_button (GtkWidget *hbox);

GtkWidget *back_button (GtkWidget *hbox);

GtkWidget *get_open_dialog (void);

void set_open_dialog (GtkWidget *w);

void edit_dialog (const char *diagtxt, const char *infotxt, const char *deftext, 
		  void (*okfunc)(), void *okptr,
		  guint hlpcode, guint varclick, int blocking);

const gchar *edit_dialog_get_text (dialog_t *dlg);

gchar *edit_dialog_special_get_text (dialog_t *dlg);

int edit_dialog_get_action (const dialog_t *dlg);

gretlopt edit_dialog_get_opt (const dialog_t *dlg);

gpointer edit_dialog_get_data (dialog_t *dlg);

void close_dialog (dialog_t *dlg);

#ifdef OLD_GTK
GtkWidget *standard_button (int code);
#endif

#endif /* DLGUTILS_H */
