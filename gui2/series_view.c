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
#include "textutil.h"
#include "dlgutils.h"
#include "clipboard.h"

typedef struct {
    char label[OBSLEN];
    double val;
} data_point_t;    

typedef struct {
    int varnum;
    int npoints;
    int digits;
    char format;    
    GtkWidget *digit_spin;
    data_point_t *points;
} series_view_t;

static void series_view_format_dialog (windata_t *vwin, guint a, GtkWidget *w);
static void series_view_sort (windata_t *vwin, guint a, GtkWidget *w);
static void series_view_graph (windata_t *vwin, guint a, GtkWidget *w);
static void scalar_to_clipboard (windata_t *vwin, guint a, GtkWidget *w);

#ifndef OLD_GTK
GtkItemFactoryEntry series_view_items[] = {
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, 
      "<StockItem>", GTK_STOCK_COPY },
    { N_("/Edit/Copy _all"), "", text_copy, COPY_TEXT, "<StockItem>", GTK_STOCK_COPY },
    { N_("/Edit/_Format..."), NULL, series_view_format_dialog, 0, NULL, GNULL },
    { N_("/_Series"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Series/_Sort"), NULL, series_view_sort, 0, NULL, GNULL },
    { N_("/Series/_Graph"), NULL, series_view_graph, 0, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

GtkItemFactoryEntry scalar_view_items[] = {
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Format..."), NULL, series_view_format_dialog, 0, NULL, GNULL },
    { N_("/Edit/_Copy value"), NULL, scalar_to_clipboard, 0, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};
#else
GtkItemFactoryEntry series_view_items[] = {
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { N_("/Edit/_Format..."), NULL, series_view_format_dialog, 0, NULL },
    { N_("/_Series"), NULL, NULL, 0, "<Branch>" },
    { N_("/Series/_Sort"), NULL, series_view_sort, 0, NULL },
    { N_("/Series/_Graph"), NULL, series_view_graph, 0, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

GtkItemFactoryEntry scalar_view_items[] = {
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Format..."), NULL, series_view_format_dialog, 0, NULL },
    { N_("/Edit/_Copy value"), NULL, scalar_to_clipboard, 0, NULL },
    { NULL, NULL, NULL, 0, NULL }
};
#endif


GtkItemFactoryEntry *get_series_view_menu_items (int code)
{
    if (code == VIEW_SERIES) return series_view_items;
    else return scalar_view_items;
}

void free_series_view (gpointer p)
{
    series_view_t *sview = (series_view_t *) p;

    if (sview == NULL) return;

    if (sview->points != NULL) {
	free(sview->points);
    }

    free(sview);
}

static int series_view_allocate (series_view_t *sview)
{
    if (sview->npoints != 0) {
	/* already allocated */
	return 0;
    }

    else if (!datainfo->vector[sview->varnum]) {
	sview->npoints = 1;
	return 0;
    }

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

static void series_view_print (windata_t *vwin)
{
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;
#endif
    PRN *prn;
    series_view_t *sview = (series_view_t *) vwin->data;
    int t;

    if (bufopen(&prn)) return;
    
    /* print formatted data to buffer */
    if (datainfo->vector[sview->varnum]) {
	pprintf(prn, "\n     Obs ");
	pprintf(prn, "%13s\n\n", datainfo->varname[sview->varnum]);
	for (t=0; t<sview->npoints; t++) {
	    if (na(sview->points[t].val)) {
		pprintf(prn, "%*s\n", OBSLEN - 1, sview->points[t].label);
	    } else if (sview->format == 'G') {
		pprintf(prn, "%*s %#13.*g\n", OBSLEN - 1, 
			sview->points[t].label,
			sview->digits, sview->points[t].val);
	    } else {
		pprintf(prn, "%*s %13.*f\n", OBSLEN - 1, 
			sview->points[t].label,
			sview->digits, sview->points[t].val);
	    }
	}
    } else {
	if (sview->format == 'G') {
	    pprintf(prn, "\n%*s = %#13.*g", OBSLEN - 1,
		    datainfo->varname[sview->varnum], 
		    sview->digits, Z[sview->varnum][0]);
	} else {
	    pprintf(prn, "\n%*s = %13.*fg", OBSLEN - 1,
		    datainfo->varname[sview->varnum], 
		    sview->digits, Z[sview->varnum][0]);
	}
    }

    /* clear existing text buffer and insert data */
#ifndef OLD_GTK
    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gtk_text_buffer_set_text(tbuf, prn->buf, -1);
#else
    gtk_text_freeze(GTK_TEXT(vwin->w));
    gtk_editable_delete_text(GTK_EDITABLE(vwin->w), 0, -1);
    gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
		    NULL, NULL, prn->buf, 
		    strlen(prn->buf));
    gtk_text_thaw(GTK_TEXT(vwin->w));
#endif

    gretl_print_destroy(prn);
}

static int compare_points (const void *a, const void *b)
{
    const data_point_t *pa = (const data_point_t *) a;
    const data_point_t *pb = (const data_point_t *) b;
     
    return (pa->val > pb->val) - (pa->val < pb->val);
}

static void series_view_sort (windata_t *vwin, guint action, GtkWidget *w)
{
    int err;
    series_view_t *sview = (series_view_t *) vwin->data;
    
    err = series_view_allocate(sview);

    if (!err) {
	/* sort the data */
	qsort((void *) sview->points, (size_t) sview->npoints, 
	      sizeof sview->points[0], compare_points);

	/* print sorted data to buffer */
	series_view_print(vwin);
    }
}

static void 
series_view_graph (windata_t *vwin, guint action, GtkWidget *w)
{
    series_view_t *sview = (series_view_t *) vwin->data;

    if (dataset_is_time_series(datainfo)) {
	do_graph_var(sview->varnum);
    } else {
	do_boxplot_var(sview->varnum);
    }
}

static void 
scalar_to_clipboard (windata_t *vwin, guint action, GtkWidget *w)
{
    series_view_t *sview = (series_view_t *) vwin->data;
    double val;
    gchar *buf;

    val = Z[sview->varnum][0];

    if (sview->format == 'G') {
	buf = g_strdup_printf("%#.*g", sview->digits, val);
    } else {
	buf = g_strdup_printf("%.*fg", sview->digits, val);
    }

#ifdef G_OS_WIN32
    win_buf_to_clipboard(buf);
#else
    buf_to_clipboard(buf);
#endif

    g_free(buf);
}

void series_view_connect (windata_t *vwin, int varnum)
{
    series_view_t *sview;

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
#ifndef OLD_GTK
        i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
#else
        i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
#endif
        sview->format = i;
    }
}

static void 
series_view_format_dialog (windata_t *vwin, guint action, GtkWidget *src)
{
    GtkWidget *w, *tmp, *label;
    GtkWidget *vbox, *hbox;
    GtkObject *adj;
    GSList *group;
    series_view_t *sview = (series_view_t *) vwin->data;

    if (series_view_allocate(sview)) return;

    w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(w), _("gretl: data format"));
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(w), "destroy",  
		     G_CALLBACK(gtk_main_quit), NULL);
#else
    gtk_signal_connect(GTK_OBJECT(w), "destroy",  
		     GTK_SIGNAL_FUNC(gtk_main_quit), NULL);
#endif

    vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 5);

    label = gtk_label_new(_("Select data format"));
    gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, FALSE, 5);

    hbox = gtk_hbox_new (FALSE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, TRUE, TRUE, 5);

    /* spinner for number of digits */
    tmp = gtk_label_new(_("Show"));
    adj = gtk_adjustment_new(sview->digits, 1, 10, 1, 1, 1);
    sview->digit_spin = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
#ifndef OLD_GTK
    g_signal_connect (adj, "value_changed",
		      G_CALLBACK (series_view_get_figures), sview);
#else
    gtk_signal_connect (adj, "value_changed",
			GTK_SIGNAL_FUNC (series_view_get_figures), sview);
#endif
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), sview->digit_spin, FALSE, FALSE, 5);

    /* select decimal places versus significant figures */
    tmp = gtk_radio_button_new_with_label (NULL, _("significant figures"));
    gtk_box_pack_start (GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    if (sview->format == 'G')
        gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (tmp), TRUE);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", 
                      GINT_TO_POINTER('G'));
#else
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(set_series_float_format), sview);
    gtk_object_set_data(GTK_OBJECT(tmp), "action", 
			GINT_TO_POINTER('G'));
#endif

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (tmp));
    tmp = gtk_radio_button_new_with_label(group, _("decimal places"));
    gtk_box_pack_start (GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    if (sview->format == 'f')
        gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (tmp), TRUE);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", 
                      GINT_TO_POINTER('f')); 
#else
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(set_series_float_format), sview);
    gtk_object_set_data(GTK_OBJECT(tmp), "action", 
			GINT_TO_POINTER('f'));    
#endif   

    /* control buttons */
    hbox = gtk_hbox_new (TRUE, 5);
    tmp = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, 0);
#ifndef OLD_GTK
    g_signal_connect_swapped (G_OBJECT (tmp), "clicked", 
			      G_CALLBACK (gtk_widget_destroy), 
			      G_OBJECT (w));
#else
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC (delete_widget), 
			GTK_OBJECT (w));
#endif

    tmp = standard_button(GTK_STOCK_CANCEL);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, 0);
#ifndef OLD_GTK
    g_signal_connect (G_OBJECT (tmp), "clicked", 
		      G_CALLBACK (series_view_format_cancel), sview);
    g_signal_connect_swapped (G_OBJECT (tmp), "clicked", 
			      G_CALLBACK (gtk_widget_destroy), 
			      G_OBJECT (w));
#else
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC (series_view_format_cancel), sview);
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC (delete_widget), 
			GTK_OBJECT (w));
#endif

    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    gtk_container_add(GTK_CONTAINER(w), vbox);

    gtk_widget_show_all(w);

#ifdef OLD_GTK
    gtk_window_set_transient_for(GTK_WINDOW(w),
				 GTK_WINDOW(vwin->dialog));
#endif

    gtk_window_set_modal(GTK_WINDOW(w), TRUE);

    gtk_main(); /* block */

    if (sview->digits > 0) {
	series_view_print(vwin);
    } else { /* canceled */
	sview->digits = 6;
    }
}

