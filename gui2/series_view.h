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

#ifndef SERIES_VIEW_H
#define SERIES_VIEW_H

#include "gretl.h"

typedef struct multi_series_view_t multi_series_view;

void free_series_view (gpointer p);

void free_multi_series_view (gpointer p);

void series_view_connect (windata_t *vwin, int varnum);

multi_series_view *multi_series_view_new (int *list);

void series_view_graph (GtkWidget *w, windata_t *vwin);

void series_view_format_dialog (GtkWidget *src, windata_t *vwin);

void series_view_sort (GtkWidget *w, windata_t *vwin);

void series_view_sort_by (GtkWidget *w, windata_t *vwin);

int *series_view_get_list (windata_t *vwin);

int series_view_is_sorted (windata_t *vwin);

PRN *vwin_print_sorted_as_csv (windata_t *vwin);

void scalar_to_clipboard (windata_t *vwin);

#endif

