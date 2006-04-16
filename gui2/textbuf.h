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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
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

gchar *textview_get_text (GtkWidget *view);

int textview_insert_text (GtkWidget *view, const gchar *text);

int viewer_char_count (windata_t *vwin);

void text_paste (windata_t *vwin, guint u, GtkWidget *widget);

void text_undo (windata_t *vwin, guint u, GtkWidget *widget);

void text_buffer_insert_colorized_buffer (GtkTextBuffer *tbuf, PRN *prn);

void text_buffer_insert_file (GtkTextBuffer *tbuf, const char *fname, 
			      int role);

GtkWidget *create_text (GtkWidget *dlg, GtkTextBuffer **buf, 
			int hsize, int vsize, gboolean editable);

void text_table_setup (GtkWidget *vbox, GtkWidget *w);

void set_help_topic_buffer (windata_t *hwin, int hcode, int pos, int en);

gboolean help_popup_handler (GtkWidget *w, GdkEventButton *event, 
			     gpointer p);

#ifdef USE_GTKSOURCEVIEW

void create_source (windata_t *vwin, GtkSourceBuffer **buf, 
		    int hsize, int vsize, gboolean editable);

void source_buffer_insert_file (GtkSourceBuffer *sbuf, const char *filename,
				int role);
#else

void correct_line_color (windata_t *vwin);

#endif

#endif /* TEXTBUF_H */
