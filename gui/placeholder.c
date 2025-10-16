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

#include "gretl.h"
#include "placeholder.h"

#if GTK_MAJOR_VERSION == 2

/* For gtk2, emulate gtk3's place-holder text apparatus: make the dummy
   text gray, and make it disappear when the GtkEntry receives focus.
*/

static gboolean undo_placeholder (GtkWidget *w, gpointer p)
{
    if (widget_get_int(w, PLACEHOLDER)) {
	if (p == NULL) {
	    /* just getting focus, not inserting text */
	    gtk_entry_set_text(GTK_ENTRY(w), "");
	}
	gtk_widget_modify_text(w, GTK_STATE_NORMAL, NULL);
	widget_set_int(w, PLACEHOLDER, 0);
    }

    return FALSE;
}

static void make_text_gray (GtkWidget *entry)
{
    GdkColor gray = {
	0, 32767, 32767, 32767
    };

    gtk_widget_modify_text(entry, GTK_STATE_NORMAL, &gray);
}

#endif /* GTK_MAJOR_VERSION == 2 */

void set_placeholder_text (GtkWidget *w, const char *s)
{
#if GTK_MAJOR_VERSION > 2
    gtk_entry_set_placeholder_text(GTK_ENTRY(w), s);
#else
    gtk_entry_set_text(GTK_ENTRY(w), s);
    make_text_gray(w);
    widget_set_int(w, PLACEHOLDER, 1);
    if (!widget_get_int(w, "signal-set")) {
	g_signal_connect(G_OBJECT(w), "grab-focus",
			 G_CALLBACK(undo_placeholder), NULL);
	g_signal_connect(G_OBJECT(w), "changed",
			 G_CALLBACK(undo_placeholder), w);
	widget_set_int(w, "signal-set", 1);
    }
#endif
}

/* equivalent to entry_box_get_trimmed_text() except that it handles the
   case where we've "faked" placeholder text for use with gtk2: such
   text is skipped.
*/

gchar *entry_box_get_real_text (GtkWidget *w)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    gchar *ret = NULL;

#if GTK_MAJOR_VERSION == 2
    if (widget_get_int(w, PLACEHOLDER)) {
        ; /* no real text available */
    }
#else
    if (s != NULL) {
	while (isspace(*s)) s++;
	if (*s != '\0') {
	    ret = g_strstrip(g_strdup(s));
	}
    }
#endif

    return ret;
}
