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

/* session.h for gretl */

#ifndef SESSION_H
#define SESSION_H

#ifndef G_OS_WIN32
# include <gtkextra/gtkiconlist.h>
#else
# include "gtkiconlist.h"
#endif 

typedef struct {
    gchar *name;
    gint sort;
    gpointer data;
    GtkIconListItem *icon;
} gui_obj;

void session_state (gboolean s);

void session_close_state (gboolean s);

void add_last_graph (gpointer data, guint code, GtkWidget *w);

void remember_model (gpointer data, guint close, GtkWidget *widget);

void session_init (void);

void do_open_session (GtkWidget *w, gpointer data);

void close_session (void);

void free_session (void);

int saved_objects (char *fname);

int parse_savefile (char *fname, SESSION *psession, session_t *rebuild);

int recreate_session (char *fname, SESSION *psession, session_t *rebuild);

void view_session (void);

#endif /* SESSION_H */
