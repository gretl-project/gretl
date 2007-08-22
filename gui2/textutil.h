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

void window_tex_callback (GtkWidget *w, windata_t *vwin);

void model_tex_view (gpointer data, guint fmt, GtkWidget *w);

void model_tex_save (gpointer data, guint fmt, GtkWidget *w);

void var_tex_callback (gpointer data, guint opt, GtkWidget *w);

void window_copy (gpointer data, guint how, GtkWidget *w);

void window_save (windata_t *vwin, guint fmt);

void text_replace (windata_t *vwin, guint u, GtkWidget *w);

#if defined(G_OS_WIN32) || defined (USE_GNOME)
void window_print (windata_t *vwin, guint u, GtkWidget *w);
#endif

void system_print_buf (const gchar *buf, FILE *fp);

char *dosify_buffer (const char *buf, int format);

#endif /* TEXTUTIL_H */
