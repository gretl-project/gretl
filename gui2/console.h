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

#ifndef CONSOLE_H
#define CONSOLE_H

void show_gretl_console (void);

gint console_key_handler (GtkWidget *w, GdkEventKey *key, 
			      gpointer p);

gint console_mouse_handler (GtkWidget *w, GdkEventButton *event,
			    gpointer p);

void console_record_sample (const DATAINFO *pdinfo);

#if 0
gint console_click_handler (GtkWidget *w, GdkEventButton *event,
			    gpointer p);
#endif


#endif /* CONSOLE_H */
