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

#include "objstack.h"

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

int save_session (char *fname);

int session_is_saved (void);

void set_session_saved (int val);

void session_menu_state (gboolean s);

const char *get_session_dirname (void);

int real_add_text_to_session (PRN *prn, const char *tname);

int add_graph_to_session (char *fname, char *fullname);

void add_boxplot_to_session (void);

int cli_add_graph_to_session (const char *fname, const char *gname,
			      GretlObjType type);

void model_add_as_icon (gpointer p, guint type, GtkWidget *w);

void model_add_as_icon_and_close (gpointer p, guint type, GtkWidget *w);

int maybe_add_model_to_session (void *ptr, GretlObjType type);

void session_model_callback (void *ptr, int action);

void *get_session_object_by_name (const char *name, GretlObjType *type);

int session_matrix_destroy_by_name (const char *name, PRN *prn);

void delete_text_from_session (void *p);

void display_saved_text (void *p);

int session_changed (int set);

void session_init (void);

void do_open_session (void);

void gui_clear_dataset (void);

void verify_clear_data (void);

void close_session (ExecState *s, double ***pZ, DATAINFO **ppdinfo);

void free_session (void);

int highest_numbered_variable_in_session (void);

int is_session_file (const char *fname);

int is_session_model (void *p);

void view_session (void);

void save_session_callback (GtkWidget *w, guint i, gpointer data);

int session_file_is_open (void);

int clear_or_save_model (MODEL **ppmod, DATAINFO *pdinfo, int rebuild);

void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w);

void save_plot_commands_callback (GtkWidget *w, gpointer p);

void disable_graph_page (void);

void display_session_graph_by_data (void *p);

void view_matrix_properties (const gretl_matrix *m, const char *name);

#endif /* SESSION_H */
