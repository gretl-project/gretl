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

#include <gtk/gtk.h>
#include "focus.h"

/* This module added 2021-06-17. The point of wrapping gtk_dialog_new()
   and gtk_window_new() is to inject common code that supports a
   record of the currently focused window. Then when functions such
   as errbox, infobox and warnbox (or variants thereof) are called
   we can arrange for the resulting message dialogs to be parented
   by the window from which they were spawned. Hopefully this will
   improve the placement of the transient windows.

   The code that utilizes get_focus_window() is to be found in
   dialogs.c, in the msgbox() function.
*/

#define FOCUS_DEBUG 0

/* the window that currently has focus */
static GtkWidget *focus_window;

GtkWidget *get_focus_window (void)
{
    return focus_window;
}

/* callback for when a window receives focus */

static gboolean focus_in_cb (GtkWidget *w, GdkEvent *event,
			     gpointer p)
{
#if FOCUS_DEBUG
    fprintf(stderr, "focus_in: %p\n", (void *) w);
#endif
    focus_window = w;
    return FALSE;
}

/* callback for when a window loses focus */

static gboolean focus_out_cb (GtkWidget *w, GdkEvent *event,
			      gpointer p)
{
#if FOCUS_DEBUG
    fprintf(stderr, "focus_out: %p\n", (void *) w);
#endif
    focus_window = NULL;
    return FALSE;
}

/* callback for window destroyed */

static void destroy_cb (GtkWidget *w, gpointer p)
{
#if FOCUS_DEBUG
    fprintf(stderr, "focus destroy: %p\n", (void *) w);
#endif
    if (focus_window == w) {
	focus_window = NULL;
    }
}

static GtkWidget *real_gretl_gtk_object (int window)
{
    GtkWidget *ret;

    if (window) {
	ret = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    } else {
	ret = gtk_dialog_new();
	/* redundant? */
	gtk_widget_add_events(ret, GDK_FOCUS_CHANGE_MASK);
    }

    /* attach our focus-related callbacks */
    g_signal_connect(G_OBJECT(ret), "focus-in-event",
		     G_CALLBACK(focus_in_cb), NULL);
    g_signal_connect(G_OBJECT(ret), "focus-out-event",
		     G_CALLBACK(focus_out_cb), NULL);
    g_signal_connect(G_OBJECT(ret), "destroy",
		     G_CALLBACK(destroy_cb), NULL);

    return ret;
}

/* wrapper for gtk_dialog_new() */

GtkWidget *gretl_gtk_dialog (void)
{
    return real_gretl_gtk_object(0);
}

/* wrapper for gtk_window_new (toplevel) */

GtkWidget *gretl_gtk_window (void)
{
    return real_gretl_gtk_object(1);
}
