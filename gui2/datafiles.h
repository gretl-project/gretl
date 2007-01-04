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

/* datafiles.h for gretl */

#ifndef DATAFILES_H
#define DATAFILES_H

void browser_open_data (GtkWidget *w, gpointer data);

void browser_open_ps (GtkWidget *w, gpointer data);

void browser_load_func (GtkWidget *w, gpointer data);

void browser_edit_func (GtkWidget *w, gpointer data);

void browser_call_func (GtkWidget *w, gpointer data);

windata_t *gui_show_function_info (const char *fname, int role);

void destroy_file_collections (void);

void display_files (gpointer p, guint code, GtkWidget *w);

gint populate_filelist (windata_t *fdata, gpointer p);

char *strip_extension (char *s);

#endif /* DATAFILES_H */
