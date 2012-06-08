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

#ifndef SSHEET_H
#define SSHEET_H

typedef enum {
    SHEET_EDIT_VARLIST,
    SHEET_EDIT_DATASET,
    SHEET_NEW_DATASET,
    SHEET_EDIT_MATRIX,
    SHEET_EDIT_SCALARS
} SheetCmd;

void show_spreadsheet (SheetCmd c);

void show_spreadsheet_for_series (int varnum);

void edit_scalars (void);

void sync_scalars_window (void);

void gui_new_matrix (GtkWidget *parent);

void edit_user_matrix_by_name (const char *name, GtkWidget *parent);

void set_dataset_locked (gboolean s);

int dataset_locked (void);

gint get_string_width (const gchar *str);

#endif /* SSHEET_H */
