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

#ifndef TEXTBUF_H
#define TEXTBUF_H

#define textview_use_highlighting(r) (r != EDIT_X12A &&			\
				      (vwin_editing_script(r) ||	\
				       r == VIEW_SCRIPT ||		\
				       r == VIEW_PKG_CODE ||		\
				       r == EDIT_PKG_CODE ||		\
				       r == EDIT_PKG_SAMPLE ||		\
				       r == VIEW_PKG_SAMPLE ||		\
				       r == VIEW_LOG ||			\
				       r == CONSOLE))

#define lbracket(k) (k == GDK_parenleft || \
		     k == GDK_bracketleft || \
		     k == GDK_braceleft)

extern int tabwidth;
extern int smarttab;
extern int script_line_numbers;
extern int script_auto_bracket;

void text_set_cursor (GtkWidget *w, GdkCursorType cspec);

void cursor_to_top (windata_t *vwin);

void cursor_to_end (windata_t *vwin);

void cursor_to_mark (windata_t *vwin, GtkTextMark *mark);

void scroll_to_foot (windata_t *vwin);

void scroll_to_line (windata_t *vwin, int line);

int textbuf_get_n_lines (windata_t *vwin);

gint get_char_width (GtkWidget *widget);

gint get_char_height (GtkWidget *widget);

int get_screen_height (void);

GtkTextBuffer *gretl_text_buf_new (void);

void gretl_viewer_set_formatted_buffer (windata_t *vwin, const char *buf);

void gretl_viewer_insert_formatted_buffer (windata_t *vwin, const char *buf);

gchar *textview_get_text (GtkWidget *view);

gchar *textview_get_hansl (GtkWidget *view, int save);

gchar *textview_get_trimmed_text (GtkWidget *view);

gchar *textview_get_normalized_line (GtkWidget *view);

gchar *textview_get_wrapped_text (GtkWidget *view);

gchar *textview_get_selection_or_all (GtkWidget *view,
				      gboolean *selection);

int textview_set_text (GtkWidget *view, const gchar *text);

void textview_append_text (GtkWidget *view, const char *text);

void textview_insert_text (GtkWidget *view, const char *text);

void textview_clear_text (GtkWidget *view);

int textview_set_text_selected (GtkWidget *view, const gchar *text);

int textview_set_cursor_at_line (GtkWidget *view, int line);

int viewer_char_count (windata_t *vwin);

void text_paste (GtkWidget *w, windata_t *vwin);

void text_redo (GtkWidget *w, windata_t *vwin);

void text_undo (GtkWidget *w, windata_t *vwin);

int text_can_undo (windata_t *vwin);

void text_larger (GtkWidget *w, gpointer data);

void text_smaller (GtkWidget *w, gpointer data);

void revise_gretl_mono_tags (void);

void textview_set_text_colorized (GtkWidget *view, const char *buf);

void textview_append_text_colorized (GtkWidget *view, const char *buf,
				     int trim);

void textview_set_text_report (GtkWidget *view, const char *buf);

void textview_set_text_dbsearch (windata_t *vwin, const char *buf);

void textview_insert_file (windata_t *vwin, const char *fname);

void create_text (windata_t *vwin, int hsize, int vsize, 
		  int nlines, gboolean editable);

void text_set_word_wrap (GtkWidget *w, gboolean wrap);

void text_table_setup (GtkWidget *vbox, GtkWidget *w);

int set_help_topic_buffer (windata_t *hwin, int pos);

gboolean help_popup_handler (GtkWidget *w, GdkEventButton *event, 
			     gpointer p);

void create_source (windata_t *vwin, int hsize, int vsize, 
		    gboolean editable);

GtkWidget *create_sample_source (const char *style);

void set_style_for_textview (GtkWidget *text, const char *id);

void create_console (windata_t *vwin, int hsize, int vsize);

void sourceview_insert_file (windata_t *vwin, const char *fname);

void sourceview_insert_buffer (windata_t *vwin, const char *buf);

void sourceview_print (windata_t *vwin);

const char **get_sourceview_style_ids (int *n);

const char **get_graph_theme_ids (int *n);

void update_script_editor_options (windata_t *vwin);

void viewer_split_pane (windata_t *vwin, int vertical);

void viewer_close_pane (windata_t *vwin);

void textview_add_processing_message (GtkWidget *view);

void textview_delete_processing_message (GtkWidget *view);

void textview_format_paragraph (GtkWidget *view);

void indent_hansl (GtkWidget *w, windata_t *vwin);

gboolean script_bracket_handler (windata_t *vwin, guint keyval);

#ifdef HAVE_GTKSV_COMPLETION

int maybe_try_completion (windata_t *vwin);

#endif

#endif /* TEXTBUF_H */
