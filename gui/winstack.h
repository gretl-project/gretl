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

#ifndef WINSTACK_H_
#define WINSTACK_H_

void window_list_add (GtkWidget *w, int role);

void window_list_popup (GtkWidget *src, GdkEvent *event,
			gpointer data);

void window_list_revise_label (GtkWidget *targ,
			       const char *label);

void plot_window_set_label (GtkWidget *w);

void vwin_record_toolbar_popup (windata_t *vwin,
				GtkWidget *menu);

void vwin_free_toolbar_popups (windata_t *vwin);

gboolean window_list_exit_check (void);

windata_t *get_editor_for_file (const char *filename);

windata_t *get_viewer_for_file (const char *filename);

windata_t *get_browser_for_database (const char *filename);

windata_t *get_browser_for_gretl_database (const char *filename);

windata_t *get_browser_for_role (int role, const char *path);

windata_t *get_viewer_for_data (const gpointer data);

GtkWidget *get_window_for_data (const gpointer data);

GtkWidget *get_window_for_plot (void *session_plot);

GtkWidget *get_viewer_for_plot (const char *filename);

gboolean is_shell_for_plotfile (GtkWidget *w,
				const char *fname);

gboolean package_being_edited (const char *pkgname,
			       GtkWidget **pw);

int get_script_output_number (void);

windata_t *get_unique_output_viewer (void);

int get_n_hansl_editor_windows (void);

GtkWidget *get_primary_hansl_window (void);

void maybe_close_window_for_user_var (const gpointer data,
				      GretlObjType otype);

void close_session_windows (gretlopt opt);

void cascade_session_windows (void);

void vwin_add_winlist (windata_t *vwin);

void window_add_winlist (GtkWidget *window, GtkWidget *hbox);

windata_t *vwin_new (int role, gpointer data);

int highest_numbered_variable_in_winstack (void);

GList *windowed_model_list (void);

windata_t *gretl_viewer_new (int role, const gchar *title,
			     gpointer data);

windata_t *
gretl_viewer_new_with_parent (windata_t *parent, int role,
			      const gchar *title,
			      gpointer data);

windata_t *gretl_browser_new (int role, const gchar *title);

GtkWidget *vwin_toplevel (windata_t *vwin);

void vwin_pack_toolbar (windata_t *vwin);

void vwin_reinstate_toolbar (windata_t *vwin);

void gretl_viewer_present (windata_t *vwin);

void gretl_viewer_destroy (windata_t *vwin);

void gretl_viewer_set_title (windata_t *vwin, const char *title);

void preset_viewer_flag (windata_flags f);

void vwin_record_action (windata_t *vwin, GtkAction *action);

void vwin_action_set_sensitive (windata_t *vwin,
				const char *name,
				gboolean s);

void set_main_winpos (int pos);

#endif
