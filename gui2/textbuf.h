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

#define help_role(r) (r == CLI_HELP || \
                      r == GUI_HELP || \
                      r == CLI_HELP_EN || \
                      r == GUI_HELP_EN)

void text_set_cursor (GtkWidget *w, GdkCursorType cspec);

void cursor_to_top (windata_t *vwin);

gint get_char_width (GtkWidget *widget);

GtkTextBuffer *gretl_text_buf_new (void);

gchar *textview_get_text (GtkWidget *view);

gchar *textview_get_selection_or_all (GtkWidget *view,
				      int *sel);

int textview_set_text (GtkWidget *view, const gchar *text);

int textview_set_cursor_at_line (GtkWidget *view, int line);

int viewer_char_count (windata_t *vwin);

void text_paste (windata_t *vwin, guint u, GtkWidget *widget);

void text_undo (windata_t *vwin, guint u, GtkWidget *widget);

void textview_set_text_colorized (GtkWidget *view, const char *buf);

void textview_append_text_colorized (GtkWidget *view, const char *buf);

void textview_insert_file (windata_t *vwin, const char *fname);

void textview_insert_from_tempfile (windata_t *vwin, PRN *prn);

void create_text (windata_t *vwin, int hsize, int vsize, 
		  gboolean editable);

void text_set_word_wrap (GtkWidget *w, gboolean wrap);

void text_table_setup (GtkWidget *vbox, GtkWidget *w);

void set_help_topic_buffer (windata_t *hwin, int hcode, int pos, int en);

gboolean help_popup_handler (GtkWidget *w, GdkEventButton *event, 
			     gpointer p);

gboolean 
script_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p);


void create_source (windata_t *vwin, int hsize, int vsize, 
		    gboolean editable);

void sourceview_insert_file (windata_t *vwin, const char *fname);

void sourceview_insert_buffer (windata_t *vwin, const char *buf);

#endif /* TEXTBUF_H */
