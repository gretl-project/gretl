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

/* series_view.c for gretl */

#include "gretl.h"

typedef struct {
    char label[9];
    double val;
} data_point_t;    

typedef struct {
    int varnum;
    int npoints;
    int digits;
    GtkWidget *digit_spin;
    char format;
    data_point_t *points;
} series_view_t;

static void series_view_format_dialog (GtkWidget *w, 
				       windata_t *vwin);

void free_series_view (gpointer p)
{
    series_view_t *sview = (series_view_t *) p;

    if (sview == NULL || sview->points == NULL) return;

    free(sview->points);
}

static int series_view_allocate (series_view_t *sview)
{
    if (sview->npoints != 0) return 0; /* already allocated */
    else {
	int t, tp, T = datainfo->t2 - datainfo->t1 + 1;
	int v = sview->varnum;

	/* allocate storage */
	sview->points = malloc(T * sizeof *sview->points);
	if (sview->points == NULL) {
	    return 1;
	} 
	sview->npoints = T;

	/* populate from data set */
	for (t=datainfo->t1; t<=datainfo->t2; t++) {
	    tp = t - datainfo->t1; 
	    sview->points[tp].val = Z[v][t];
	    if (datainfo->markers) {
		strcpy(sview->points[tp].label, datainfo->S[t]);
	    } else {
		ntodate(sview->points[tp].label, t, datainfo);
	    }
	}
    }

    return 0;
}

static int compare_points (const void *a, const void *b)
{
    const data_point_t *pa = (const data_point_t *) a;
    const data_point_t *pb = (const data_point_t *) b;
     
    return (pa->val > pb->val);
}

static void series_view_sort (GtkWidget *w, gpointer data)
{
    int t, err;
    windata_t *vwin = (windata_t *) data;
    series_view_t *sview = (series_view_t *) vwin->data;
    PRN *prn;

    if (bufopen(&prn)) return;
    
    err = series_view_allocate(sview);

    if (!err) {
	GtkTextBuffer *tbuf;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

	/* sort the data */
	qsort(sview->points, sview->npoints, 
	      sizeof *sview->points, compare_points);

	/* print sorted data to buffer */
	pprintf(prn, "\n     Obs ");
	pprintf(prn, "%13s\n\n", datainfo->varname[sview->varnum]);
	for (t=0; t<sview->npoints; t++) {
	    pprintf(prn, "%8s %#13.6g\n", sview->points[t].label, 
		    sview->points[t].val);
	}

	/* clear existing text buffer and insert sorted data */
	gtk_text_buffer_set_text(tbuf, prn->buf, -1);
    }

    gretl_print_destroy(prn);
}

void series_view_print (windata_t *vwin)
{
    PRN *prn;
    GtkTextBuffer *tbuf;
    series_view_t *sview = (series_view_t *) vwin->data;
    int t;

    if (bufopen(&prn)) return;
    
    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

    /* print formatted data to buffer */
    pprintf(prn, "\n     Obs ");
    pprintf(prn, "%13s\n\n", datainfo->varname[sview->varnum]);
    for (t=0; t<sview->npoints; t++) {
	if (sview->format == 'G') {
	    pprintf(prn, "%8s %#13.*g\n", sview->points[t].label,
		    sview->digits, sview->points[t].val);
	} else {
	    pprintf(prn, "%8s %13.*f\n", sview->points[t].label,
		    sview->digits, sview->points[t].val);
	}
    }

    /* clear existing text buffer and insert sorted data */
    gtk_text_buffer_set_text(tbuf, prn->buf, -1);

    gretl_print_destroy(prn);
}

void build_series_view_popup (windata_t *win)
{
    if (win->popup != NULL) return;

    win->popup = gtk_menu_new();

    add_popup_item(_("Sort values"), win->popup, 
		   G_CALLBACK(series_view_sort), 
		   win);
    add_popup_item(_("Format values"), win->popup, 
		   G_CALLBACK(series_view_format_dialog), 
		   win);
}

void series_view_connect (windata_t *vwin, int varnum)
{
    series_view_t *sview;

    if (!datainfo->vector[varnum]) return;

    sview = malloc(sizeof *sview);
    if (sview == NULL) {
	vwin->data = NULL;
    } else {
	sview->varnum = varnum;
	sview->npoints = 0;
	sview->points = NULL;
	sview->digits = 6;
	sview->digit_spin = NULL;
	sview->format = 'G';
	vwin->data = sview;
    }
}

static 
void series_view_format_cancel (GtkWidget *w, series_view_t *sview)
{
    sview->digits = -1;
}

static 
void series_view_get_figures (GtkWidget *w, series_view_t *sview)
{
    sview->digits = gtk_spin_button_get_value_as_int
	(GTK_SPIN_BUTTON(sview->digit_spin));
}

static void 
set_series_float_format (GtkWidget *w, gpointer p)
{
    gint i;
    series_view_t *sview = (series_view_t *) p;

    if (GTK_TOGGLE_BUTTON(w)->active) {
        i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
        sview->format = i;
    }
}

static void series_view_format_dialog (GtkWidget *src, windata_t *vwin)
{
    GtkWidget *w, *tmp, *label;
    GtkWidget *vbox, *hbox;
    GtkObject *adj;
    GSList *group;
    series_view_t *sview = (series_view_t *) vwin->data;

    if (series_view_allocate(sview)) return;

    w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(w), _("gretl: data format"));
    g_signal_connect(G_OBJECT(w), "destroy",  
		     G_CALLBACK(gtk_main_quit), NULL);

    vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 5);

    label = gtk_label_new(_("Select data format"));
    gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, FALSE, 5);

    hbox = gtk_hbox_new (FALSE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, TRUE, TRUE, 5);

    tmp = gtk_label_new(_("figures:"));
    adj = gtk_adjustment_new(sview->digits, 1, 10, 1, 1, 1);
    sview->digit_spin = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    g_signal_connect (adj, "value_changed",
		      G_CALLBACK (series_view_get_figures), sview);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), sview->digit_spin, FALSE, FALSE, 5);

    /* decimal places versus significant figures */
    tmp = gtk_radio_button_new_with_label (NULL, _("significant figures"));
    gtk_box_pack_start (GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    if (sview->format == 'G')
        gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (tmp), TRUE);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", 
                      GINT_TO_POINTER('G'));

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (tmp));
    tmp = gtk_radio_button_new_with_label(group, _("decimal places"));
    gtk_box_pack_start (GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    if (sview->format == 'f')
        gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (tmp), TRUE);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", 
                      GINT_TO_POINTER('f'));    

    /* control buttons */
    hbox = gtk_hbox_new (TRUE, 5);
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, 0);
    g_signal_connect_swapped (G_OBJECT (tmp), "clicked", 
			      G_CALLBACK (gtk_widget_destroy), 
			      G_OBJECT (w));

    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, 0);
    g_signal_connect (G_OBJECT (tmp), "clicked", 
		      G_CALLBACK (series_view_format_cancel), sview);
    g_signal_connect_swapped (G_OBJECT (tmp), "clicked", 
			      G_CALLBACK (gtk_widget_destroy), 
			      G_OBJECT (w));

    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    gtk_container_add(GTK_CONTAINER(w), vbox);

    gtk_widget_show_all(w);

    gtk_window_set_modal(GTK_WINDOW(w), TRUE);

    gtk_main(); /* block */

    if (sview->digits > 0) {
	series_view_print(vwin);
    } else { /* canceled */
	sview->digits = 6;
    }
}

