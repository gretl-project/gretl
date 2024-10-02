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

static gboolean obs_button_input (GtkSpinButton *spin,
				  gdouble *new_val,
				  gpointer p)
{
    const gchar *obs = gtk_entry_get_text(GTK_ENTRY(spin));
    int n;

    if (g_object_get_data(G_OBJECT(spin), "newdata")) {
	n = merge_dateton(obs, (DATASET *) p);
    } else {
	n = dateton(obs, (DATASET *) p);
    }

    *new_val = n;

    return TRUE;
}

static gboolean obs_button_output (GtkSpinButton *spin, gpointer p)
{
    gpointer rset;
    gchar buf[OBSLEN];
    int n = gtk_spin_button_get_value_as_int(spin);

    ntolabel(buf, n, (DATASET *) p);

    if (strcmp(buf, gtk_entry_get_text(GTK_ENTRY(spin)))) {
	gtk_entry_set_text(GTK_ENTRY(spin), buf);
    }

    rset = g_object_get_data(G_OBJECT(spin), "rset");
    if (rset != NULL) {
	update_obs_label(NULL, rset);
    }

    return TRUE;
}

GtkWidget *obs_button_new (GtkAdjustment *adj, DATASET *dset,
			   ObsButtonRole role)
{
    GtkWidget *spin;
    int n = strlen(dset->endobs);

    spin = gtk_spin_button_new(adj, 1, 0);
    gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(spin), FALSE);
    gtk_entry_set_width_chars(GTK_ENTRY(spin), (n < 2)? 2 : n);
    gtk_spin_button_set_snap_to_ticks(GTK_SPIN_BUTTON(spin), TRUE);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(spin),
				      GTK_UPDATE_IF_VALID);
#if GTK_MAJOR_VERSION == 3
    /* remedy required for gtk3 */
    gtk_entry_set_max_width_chars(GTK_ENTRY(spin), (n < 2)? 2 : n);
#endif

    g_signal_connect(G_OBJECT(spin), "input",
		     G_CALLBACK(obs_button_input), dset);
    g_signal_connect(G_OBJECT(spin), "output",
		     G_CALLBACK(obs_button_output), dset);

    if (role) {
	g_object_set_data(G_OBJECT(spin), "role", GINT_TO_POINTER(role));
    }

    return spin;
}

GtkWidget *data_start_button (GtkAdjustment *adj, DATASET *dset)
{
    GtkWidget *spin = obs_button_new(adj, dset, 0);

    g_object_set_data(G_OBJECT(spin), "newdata", GINT_TO_POINTER(1));

    return spin;
}

int obs_button_get_value (GtkWidget *button)
{
    return gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(button));
}

const gchar *obs_button_get_string (GtkWidget *button)
{
    return gtk_entry_get_text(GTK_ENTRY(button));
}

static void alert_partner (GtkWidget *b, GtkWidget *partner)
{
    GtkSpinButton *b1, *b2;
    int b_role, t1, t2;

    b_role = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(b), "role"));

    if (b_role == OBS_BUTTON_T1) {
	b1 = GTK_SPIN_BUTTON(b);
	b2 = GTK_SPIN_BUTTON(partner);
    } else {
	b1 = GTK_SPIN_BUTTON(partner);
	b2 = GTK_SPIN_BUTTON(b);
    }

    t1 = gtk_spin_button_get_value_as_int(b1);
    t2 = gtk_spin_button_get_value_as_int(b2);

    if (t2 < t1) {
	if (b_role == OBS_BUTTON_T1) {
	    /* force t2 to adjust upward */
	    gtk_spin_button_set_value(b2, (gdouble) t1);
	} else {
	    /* force t1 to adjust downward */
	    gtk_spin_button_set_value(b1, (gdouble) t2);
	}
    }
}

void obs_button_set_partner (GtkWidget *button, GtkWidget *partner)
{
    g_signal_connect(G_OBJECT(button), "value-changed",
		     G_CALLBACK(alert_partner), partner);
}
