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

#ifndef TEXTUTIL_H
#define TEXTUTIL_H

enum {
    W_PREVIEW,
    W_COPY,
    W_SAVE
};

int prn_to_clipboard (PRN *prn, int copycode);

int multiple_formats_ok (windata_t *vwin);

void window_tex_callback (GtkWidget *w, windata_t *vwin);

void model_tex_view (GtkAction *action, gpointer data);

void model_tex_save (GtkAction *action, gpointer data);

void model_tex_copy (GtkAction *action, gpointer data);

void window_copy (windata_t *vwin, guint fmt);

void window_save (windata_t *vwin, guint fmt);

void text_replace (GtkWidget *w, windata_t *vwin);

void window_print (GtkAction *action, windata_t *vwin);

void system_print_buf (const gchar *buf, FILE *fp);

char *dosify_buffer (const char *buf, int format);

char *strip_unicode_minus (char *s);

int has_unicode_minus (const unsigned char *s);

#endif /* TEXTUTIL_H */
