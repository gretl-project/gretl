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

enum {
    ADD_OBJECT_OK,
    ADD_OBJECT_REPLACE,
    ADD_OBJECT_FAIL
};

int session_is_saved (void);

void set_session_saved (int val);

void session_menu_state (gboolean s);

int real_add_graph_to_session (const char *fname, const char *grname,
			       int code);

int real_add_text_to_session (PRN *prn, const char *tname);

void add_graph_to_session (gpointer data, guint code, GtkWidget *w);

void remember_model (gpointer data, guint close, GtkWidget *widget);

void remember_var (gpointer data, guint close, GtkWidget *widget);

void remember_sys (gpointer data, guint close, GtkWidget *widget);

int model_already_saved (const MODEL *pmod);

int try_add_model_to_session (MODEL *pmod);

int try_add_var_to_session (GRETL_VAR *var);

void *get_session_object_by_name (const char *name, int *which);

void delete_model_from_session (MODEL *pmod);

void delete_var_from_session (GRETL_VAR *var);

void delete_text_from_session (const char *tname);

void display_text_by_name (const char *tname);

int session_changed (int set);

void session_init (void);

void do_open_session (GtkWidget *w, gpointer data);

void verify_clear_data (void);

void close_session (void);

void free_session (void);

int highest_numbered_variable_in_session (void);

int saved_objects (const char *fname);

int parse_savefile (const char *fname);

int recreate_session (const char *fname);

void view_session (void);

void save_session_callback (GtkWidget *w, guint i, gpointer data);

void session_file_manager (int action, const char *fname);

int session_file_is_open (void);

int clear_or_save_model (MODEL **ppmod, DATAINFO *pdinfo, int rebuild);

void print_saved_object_specs (const char *session_base, FILE *fp);

int print_session_notes (const char *fname);

void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w);

void save_plot_commands_callback (GtkWidget *w, gpointer p);

void disable_graph_page (void);

#endif /* SESSION_H */
