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

/* session.h for gretl */

#ifndef SESSION_H
#define SESSION_H

#include "objstack.h"

typedef struct SESSION_GRAPH_ SESSION_GRAPH;

enum {
    SCHEDULE_FOR_DELETION,
    REALLY_DELETE_ALL,
    CLEAR_DELFILES
};

enum {
    ADD_OBJECT_OK,
    ADD_OBJECT_REPLACE,
    ADD_OBJECT_FAIL
};

enum {
    LOG_SAVE,
    LOG_SAVE_AS,
    LOG_OPEN,
    LOG_CLOSE,
    LOG_NULL
};

int save_session (char *fname);

int save_session_commands (char *fname);

int save_session_dataset (void);

int session_is_modified (void);

int session_is_open (void);

void set_commands_recorded (void);

int get_commands_recorded (void);

void session_menu_state (gboolean s);

int have_session_objects (void);

int widget_is_iconview (GtkWidget *w);

const char *get_session_dirname (void);

int real_add_text_to_session (PRN *prn, int pos, const char *tname);

void save_output_as_text_icon (windata_t *vwin);

int gui_add_graph_to_session (char *fname, char *fullname,
                              int type, SESSION_GRAPH **pgrf);

int cli_add_graph_to_session (const char *fname, const char *gname,
			      GretlObjType type, int display);

char *session_graph_make_path (char *path, const char *fname);

void view_plot_commands (SESSION_GRAPH *graph);

const char *last_session_graph_name (void);

void model_add_as_icon (GtkAction *action, gpointer p);

int add_model_to_session_callback (void *ptr, GretlObjType type,
				   gretlopt opt);

void session_model_callback (void *ptr, int action);

void bundle_add_as_icon (GtkAction *action, gpointer p);

void *get_session_object_by_name (const char *name, GretlObjType *type);

int session_user_var_destroy_by_name (const char *name,
				      GretlObjType type);

void delete_text_from_session (void *p);

void display_saved_text (void *p);

void mark_session_changed (void);

void session_init (void);

gboolean do_open_session (void);

gretl_bundle *open_session_as_bundle (void);

void gui_clear_dataset (void);

void verify_clear_data (void);

void close_session (gretlopt opt);

void free_session (int on_exit);

int highest_numbered_variable_in_session (void);

GList *session_model_list (void);

int is_session_model (void *p);

void view_session (void);

void maybe_view_session (void);

void maybe_sensitize_iconview (void);

void save_session_callback (GtkAction *action);

int session_file_is_open (void);

int clear_or_save_model (MODEL **ppmod, DATASET *pdinfo, int rebuild);

void disable_graph_page (void);

void display_session_graph_by_data (void *p);

gchar *session_graph_get_filename (void *p);

void view_matrix_properties (const gretl_matrix *m, const char *name);

void session_notes_callback (GtkWidget *w, gpointer p);

#endif /* SESSION_H */
