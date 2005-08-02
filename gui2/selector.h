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

#ifndef SELECTOR_H
#define SELECTOR_H

typedef struct _selector selector;

void clear_selector (void);

void delete_selection_dialog (selector *sr);

void selection_dialog (const char *title, void (*okfunc)(), guint cmdcode,
		       int preselect);

void simple_selection (const char *title, void (*okfunc)(), guint cmdcode,
		       gpointer p);

char *main_window_selection_as_string (void);

void data_save_selection_wrapper (int file_code);

int selector_code (const selector *sr);

const char *selector_list (const selector *sr);

int selector_list_hasconst (const selector *sr);

gpointer selector_get_data (const selector *sr);

gretlopt selector_get_opts (const selector *sr);

int selector_error (const selector *sr);

void maybe_clear_selector (const int *dlist);

#endif /* SELECTOR_H */
