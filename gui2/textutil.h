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

#ifndef TEXTUTIL_H
#define TEXTUTIL_H

int prn_to_clipboard (PRN *prn, int copycode);

void text_copy (gpointer data, guint how, GtkWidget *w);

void text_replace (windata_t *mydata, guint u, GtkWidget *widget);

#if defined(G_OS_WIN32) || defined (USE_GNOME)
void window_print (windata_t *vwin, guint u, GtkWidget *widget);
#endif

void system_print_buf (const gchar *buf, FILE *fp);

char *dosify_buffer (const char *buf, int format);

#endif /* TEXTUTIL_H */
