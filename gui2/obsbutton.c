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
#include "obsbutton.h"

/* in dialogs.c */
extern gboolean update_obs_label (GtkComboBox *box, gpointer data);

static gboolean obs_button_input (GtkSpinButton *spin, 
				  gdouble *new_val,
				  gpointer p)
{
    const gchar *obs = gtk_entry_get_text(GTK_ENTRY(spin));
    int n = dateton(obs, (DATAINFO *) p);

    *new_val = n;

    return TRUE;
}

static gboolean obs_button_output (GtkSpinButton *spin, gpointer p)
{
    gpointer rset;
    gchar buf[OBSLEN];
    int n;

    n = gtk_spin_button_get_value_as_int(spin);
    ntodate(buf, n, (DATAINFO *) p);

    if (strcmp(buf, gtk_entry_get_text(GTK_ENTRY(spin)))) {
	gtk_entry_set_text(GTK_ENTRY(spin), buf);
    }

    rset = g_object_get_data(G_OBJECT(spin), "rset");
    if (rset != NULL) {
	update_obs_label(NULL, rset);
    }

    return TRUE;
}

GtkWidget *obs_button_new (GtkAdjustment *adj, DATAINFO *pdinfo) 
{
    GtkWidget *spinner;
    int n = strlen(pdinfo->endobs);

    spinner = gtk_spin_button_new(adj, 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spinner), FALSE);
    gtk_entry_set_width_chars(GTK_ENTRY(spinner), (n < 2)? 2 : n);

    g_signal_connect(G_OBJECT(spinner), "input",
		     G_CALLBACK(obs_button_input), pdinfo);
    g_signal_connect(G_OBJECT(spinner), "output",
		     G_CALLBACK(obs_button_output), pdinfo);

    return spinner;
}

int obs_button_get_value (GtkWidget *button)
{
    return gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(button));
}
