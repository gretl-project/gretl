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

#ifndef TABWIN_H_
#define TABWIN_H_

typedef struct tabwin_t_  tabwin_t;

windata_t *viewer_tab_new (int role, const char *info,
			   gpointer data);

void tabwin_register_toolbar (windata_t *vwin);

void tabwin_tab_set_title (windata_t *vwin, const char *title);

void tabwin_tab_set_status (windata_t *vwin);

void show_tabbed_viewer (windata_t *vwin);

void maybe_destroy_tabwin (windata_t *vwin);

void tabwin_navigate (windata_t *vwin, guint key);

void undock_tabbed_viewer (GtkWidget *w, windata_t *vwin);

gboolean window_is_undockable (windata_t *vwin);

void add_undock_popup_item (GtkWidget *menu, windata_t *vwin);

gboolean window_is_dockable (windata_t *vwin);

void add_dock_popup_item (GtkWidget *menu, windata_t *vwin);

gboolean tabwin_exit_check (GtkWidget *w);

void tabwin_tab_destroy (windata_t *vwin);

windata_t *tabwin_get_editor_for_file (const char *filename,
				       GtkWidget *w);

void tabwin_tab_present (windata_t *vwin);

void tabwin_close_models_viewer (GtkWidget *w);

void tabwin_register_dialog (GtkWidget *w, gpointer p);

void script_editor_show_new_open (windata_t *vwin, gboolean show);

int viewer_n_siblings (windata_t *vwin);

int highest_numbered_var_in_tabwin (tabwin_t *tabwin, 
				    const DATASET *dset);

windata_t *window_get_active_vwin (GtkWidget *window);

#endif /* TABWIN_H_ */
