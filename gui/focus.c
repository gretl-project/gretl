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

static GtkWidget *focus_window;

GtkWidget *get_focus_window (void)
{
    return focus_window;
}

static gboolean focus_in_cb (GtkWidget *w, GdkEvent *event,
			     gpointer p)
{
    fprintf(stderr, "focus_in: %p\n", (void *) w);
    return FALSE;
}

static gboolean focus_out_cb (GtkWidget *w, GdkEvent *event,
			      gpointer p)
{
    fprintf(stderr, "focus_out: %p\n", (void *) w);
    return FALSE;
}

static GtkWidget *real_gretl_gtk_object (int window)
{
    GtkWidget *ret;

    if (window) {
	ret = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    } else {
	ret = gtk_dialog_new();
    }

    gtk_widget_set_events(ret, GDK_FOCUS_CHANGE_MASK);
    g_signal_connect(G_OBJECT(ret), "focus-in-event",
		     G_CALLBACK(focus_in_cb), NULL);
    g_signal_connect(G_OBJECT(ret), "focus-out-event",
		     G_CALLBACK(focus_out_cb), NULL);

    return ret;
}

GtkWidget *gretl_gtk_dialog (void)
{
    return real_gretl_gtk_object(0);
}

GtkWidget *gretl_gtk_window (void)
{
    return real_gretl_gtk_object(1);
}
