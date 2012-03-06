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

windata_t *editor_tab_new (const char *filename);

void tabwin_register_toolbar (windata_t *vwin);

void tabwin_set_tab_title (windata_t *vwin, gchar *fname);

void show_tabbed_viewer (GtkWidget *vmain);

void maybe_destroy_tabwin (windata_t *vwin);

#endif /* TABWIN_H_ */
