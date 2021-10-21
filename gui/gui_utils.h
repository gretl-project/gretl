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

#ifndef GUI_UTILS_H
#define GUI_UTILS_H

#define vwin_editing_script(r) (r == EDIT_HANSL ||  \
				r == EDIT_GP ||	    \
				r == EDIT_R ||	    \
				r == EDIT_OX ||     \
                                r == EDIT_OCTAVE || \
				r == EDIT_PYTHON || \
				r == EDIT_STATA ||  \
				r == EDIT_JULIA ||  \
				r == EDIT_DYNARE || \
				r == EDIT_LPSOLVE || \
				r == EDIT_X12A ||   \
				r == EDIT_SPEC)

#define vwin_editing_buffer(r) (r == EDIT_HEADER || \
	                        r == EDIT_NOTES)

#define vwin_content_changed(v) (v->flags & VWIN_CONTENT_CHANGED)

#ifdef G_OS_WIN33
#define IS_SLASH(c) (c == '\\' || c == '/')
#else
#define IS_SLASH(c) (c == '/')
#endif

enum { ADD_COMMA, DROP_COMMA };

int vwin_is_editing (windata_t *vwin);

int copyfile (const char *src, const char *dest);

FILE *gretl_tempfile_open (char *fname);

void delete_widget (GtkWidget *widget, gpointer data);

void set_wait_cursor (GdkWindow **pcwin);

void unset_wait_cursor (GdkWindow *cwin);

void *mymalloc (size_t size); 

void *myrealloc (void *ptr, size_t size);

void nomem (void);

void mark_dataset_as_modified (void);

void register_data (int flags);

void register_startup_data (const char *fname);

gboolean do_open_data (windata_t *vwin, int code);

gboolean verify_open_data (windata_t *vwin, int code);

gboolean verify_open_session (void);

void free_windata (GtkWidget *w, gpointer data);

void mark_vwin_content_saved (windata_t *vwin);

void mark_vwin_content_changed (windata_t *vwin);

void vwin_set_filename (windata_t *vwin, const char *fname);

int vwin_subselection_present (windata_t *vwin);

int latex_is_ok (void);

gint query_save_text (GtkWidget *w, GdkEvent *event, windata_t *vwin);

void vwin_save_callback (GtkWidget *w, windata_t *vwin);

gboolean vwin_copy_callback (GtkWidget *w, windata_t *vwin);

void vwin_add_child (windata_t *parent, windata_t *child);

windata_t *view_buffer (PRN *prn, int hsize, int vsize, 
			const char *title, int role,
			gpointer data);

windata_t *view_buffer_with_parent (windata_t *parent, PRN *prn, 
				    int hsize, int vsize, 
				    const char *title, int role, 
				    gpointer data);

windata_t *view_file (const char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role);

windata_t *
view_file_with_title (const char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role, 
		      const char *given_title);

windata_t *view_formatted_text_buffer (const gchar *title, 
				       const char *buf, 
				       int hsize, int vsize,
				       int role);

windata_t *hansl_output_viewer_new (PRN *prn, int mode,
				    const char *title);

void set_reuseable_output_window (int policy, windata_t *vwin);

gchar *gretl_window_title (const char *s);

gchar *title_from_filename (const char *fname, int role, gboolean prepend);

windata_t *console_window (int hsize, int vsize);

windata_t *view_help_file (const char *filename, int role);

windata_t *edit_buffer (char **pbuf, int hsize, int vsize, 
			char *title, int role);

windata_t *vwin_first_child (windata_t *vwin);

windata_t *view_model (PRN *prn, MODEL *pmod, char *title);

windata_t *view_script (const char *filename, int editable, 
			int role);

gint catch_viewer_key (GtkWidget *w, GdkEventKey *key, 
		       windata_t *vwin);

void set_model_save_state (windata_t *vwin, gboolean s);

int gui_validate_varname_strict (const char *name, GretlType type,
				 GtkWidget *parent);

int gui_validate_varname (const char *varname, GretlType type,
			  GtkWidget *parent);

gretlopt get_tex_eqn_opt (void);

gint popup_menu_handler (GtkWidget *widget, GdkEventButton *event,
			 gpointer data);

void add_popup_item (const gchar *label, GtkWidget *menu,
		     GCallback callback, gpointer data);

gboolean text_popup_handler (GtkWidget *w, GdkEventButton *event, 
			     gpointer p);

GtkWidget *make_bundle_content_menu (windata_t *vwin);

GtkWidget *make_bundle_plot_menu (windata_t *vwin);

GtkWidget *make_bundle_save_menu (windata_t *vwin);

void *gui_get_plugin_function (const char *funcname);

int get_imported_data (char *fname, int ftype, int append);

const char *path_last_slash_const (const char *path);

char *gretl_basename (char *dest, const char *src, int addscore);

char *double_underscores (char *targ, const char *src);

gchar *double_underscores_new (const char *src);

char *adjust_fontspec_string (char *targ, const char *src,
			      int mod);

void verbose_gerror_report (GError *gerr, const char *src);

int gretl_file_get_contents (const gchar *fname, gchar **contents, gsize *size);

void start_R (const char *buf, int send_data, int interactive);

void run_foreign_script (gchar *buf, int lang, gretlopt opt);

const char *print_today (void);

#ifndef G_OS_WIN32
int browser_open (const char *url);
#endif

#endif /* GUI_UTILS_H */
