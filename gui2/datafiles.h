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

/* datafiles.h for gretl */

#ifndef DATAFILES_H
#define DATAFILES_H

void browser_open_data (GtkWidget *w, gpointer data);

void browser_open_ps (GtkWidget *w, gpointer data);

void browser_load_func (GtkWidget *w, gpointer data);

void browser_edit_func (GtkWidget *w, gpointer data);

void browser_call_func (GtkWidget *w, gpointer data);

void destroy_file_collections (void);

void show_files (GtkAction *action, gpointer p);

void display_files (int code, gpointer p);

gint populate_filelist (windata_t *fdata, gpointer p);

char *strip_extension (char *s);

windata_t *gui_show_function_info (const char *fname, int role);

void maybe_update_func_files_window (int code);

void set_funcs_dir_callback (windata_t *vwin, char *path);

windata_t *get_local_viewer (int remote_role);

#endif /* DATAFILES_H */
