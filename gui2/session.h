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

enum {
    SCHEDULE_FOR_DELETION,
    REALLY_DELETE_ALL,
    CLEAR_DELFILES
};

enum {
    SAVE_AS_IS,
    SAVE_RENAME
};

void session_menu_state (gboolean s);

int real_add_graph_to_session (const char *fname, const char *grname,
			       int code);

void add_graph_to_session (gpointer data, guint code, GtkWidget *w);

void remember_model (gpointer data, guint close, GtkWidget *widget);

int try_add_model_to_session (MODEL *pmod);

MODEL *get_session_model_by_name (const char *modname);

void delete_model_from_session (MODEL *pmod);

int session_changed (int set);

void session_init (void);

void do_open_session (GtkWidget *w, gpointer data);

void verify_clear_data (void);

void close_session (void);

void free_session (void);

int saved_objects (char *fname);

int parse_savefile (char *fname, SESSION *psession, SESSIONBUILD *rebuild);

int recreate_session (char *fname, SESSION *psession, SESSIONBUILD *rebuild);

void view_session (void);

void save_session_callback (GtkWidget *w, guint i, gpointer data);

void session_file_manager (int action, const char *fname);

int session_file_is_open (void);

#endif /* SESSION_H */
