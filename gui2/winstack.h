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

void winstack_init (void);

void winstack_destroy (void);

void winstack_add (GtkWidget *w);

void winstack_remove (GtkWidget *w);

int winstack_match_data (const gpointer p);

GtkWidget *match_window_by_data (const gpointer p);

GtkWidget *match_window_by_filename (const char *fname);

int highest_numbered_variable_in_winstack (void);

windata_t *gretl_viewer_new (int role, const gchar *title, 
			     gpointer data, int record);

windata_t *gretl_browser_new (int role, const gchar *title,
			      int record);

#endif

