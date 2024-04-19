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

#ifndef GUI_UTILS_H
#define GUI_UTILS_H

#ifdef G_OS_WIN33
#define IS_SLASH(c) (c == '\\' || c == '/')
#else
#define IS_SLASH(c) (c == '/')
#endif

int latex_is_ok (void);

void mark_dataset_as_modified (void);

void register_data (int flags);

void register_startup_data (const char *fname);

gboolean do_open_data (windata_t *vwin, int code);

gboolean verify_open_data (windata_t *vwin, int action,
			   gboolean dnd);

gboolean verify_open_session (void);

windata_t *console_window (int hsize, int vsize);

windata_t *view_model (PRN *prn, MODEL *pmod, char *title);

int gui_validate_varname_strict (const char *name, GretlType type,
				 GtkWidget *parent);

int gui_validate_varname (const char *varname, GretlType type,
			  GtkWidget *parent);

int gui_validate_varname_easy (const char *varname, GretlType type);

gretlopt get_tex_eqn_opt (void);

gint popup_menu_handler (GtkWidget *widget, GdkEventButton *event,
			 gpointer data);

void add_popup_item (const gchar *label, GtkWidget *menu,
		     GCallback callback, gpointer data);

void add_system_ui_to_vwin (windata_t *vwin);

int get_imported_data (char *fname, int ftype, int append);

GtkWidget *make_bundle_save_menu (windata_t *vwin);

#endif /* GUI_UTILS_H */
