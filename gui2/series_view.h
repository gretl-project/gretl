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

#ifndef SERIES_VIEW_H
#define SERIES_VIEW_H

#include "gretl.h"

typedef struct series_view_t series_view;

void free_series_view (gpointer p);

void free_multi_series_view (gpointer p);

void series_view_connect (windata_t *vwin, int varnum);

series_view *multi_series_view_new (const int *list);

void series_view_graph (GtkWidget *w, windata_t *vwin);

void series_view_edit (GtkWidget *w, windata_t *vwin);

void series_view_refresh (GtkWidget *w, windata_t *vwin);

void series_view_format_dialog (windata_t *vwin);

void series_view_toggle_sort (GtkWidget *w, windata_t *vwin);

void multi_series_view_sort_by (GtkWidget *w, windata_t *vwin);

int *series_view_get_list (windata_t *vwin);

int series_view_is_sorted (windata_t *vwin);

PRN *vwin_print_sorted_with_format (windata_t *vwin, PrnFormat fmt);

void scalar_to_clipboard (windata_t *vwin);

int has_sortable_data (windata_t *vwin);

int can_format_data (windata_t *vwin);

#endif

