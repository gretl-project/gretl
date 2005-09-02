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
#include "menustate.h"
#include "series_view.h"

#ifdef G_OS_WIN32
#include "gretlwin32.h"
#else
#include "clipboard.h"
#endif

typedef struct data_point_t data_point;
typedef struct multi_point_t multi_point;
typedef struct series_view_t series_view;

struct data_point_t {
    char label[OBSLEN];
    double val;
};

struct multi_point_t {
    int obsnum;
    double val;
};  

struct series_view_t {
    int varnum;
    int npoints;
    int digits;
    char format;    
    GtkWidget *digit_spin;
    data_point *points;
};

struct multi_series_view_t {
    int *list;
    int sortvar;
    int npoints;
    multi_point *points;
};

static void series_view_format_dialog (windata_t *vwin, guint a, GtkWidget *w);
static void series_view_sort (windata_t *vwin, guint a, GtkWidget *w);
static void series_view_graph (windata_t *vwin, guint a, GtkWidget *w);
static void scalar_to_clipboard (windata_t *vwin, guint a, GtkWidget *w);

#ifndef OLD_GTK

GtkItemFactoryEntry series_view_items[] = {
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy selection"), NULL, window_copy, GRETL_FORMAT_SELECTION, 
      "<StockItem>", GTK_STOCK_COPY },
    { N_("/Edit/Copy _all"), "", window_copy, GRETL_FORMAT_TXT, 
      "<StockItem>", GTK_STOCK_COPY },
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
    { N_("/Edit/_Copy selection"), NULL, window_copy, GRETL_FORMAT_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, window_copy, GRETL_FORMAT_TXT, NULL },
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
    if (code == VIEW_SERIES) {
	return series_view_items;
    } else {
	return scalar_view_items;
    } 
}

void free_series_view (gpointer p)
{
    series_view *sview = (series_view *) p;

    if (sview == NULL) return;

    if (sview->points != NULL) {
	free(sview->points);
    }

    free(sview);
}

void free_multi_series_view (gpointer p)
{
    multi_series_view *mview = (multi_series_view *) p;

    if (mview == NULL) return;

    if (mview->list != NULL) free(mview->list);
    if (mview->points != NULL) free(mview->points);

    free(mview);
}

static int series_view_allocate (series_view *sview)
{
    if (sview->npoints != 0) {
	/* already allocated */
	return 0;
    } else if (!datainfo->vector[sview->varnum]) {
	sview->npoints = 1;
	return 0;
    } else {
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

static int multi_series_view_allocate (multi_series_view *mview)
{
    if (mview->npoints != 0) {
	/* already allocated */
	return 0;
    } else {
	int T = datainfo->t2 - datainfo->t1 + 1;

	/* allocate storage */
	mview->points = mymalloc(T * sizeof *mview->points);
	if (mview->points == NULL) {
	    return 1;
	} 
	mview->npoints = T;
    }

    return 0;
}

static void mview_fill_points (multi_series_view *mview)
{
    int t, tp = 0;

    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	mview->points[tp].obsnum = t;
	mview->points[tp].val = Z[mview->sortvar][t];
	tp++;
    }
}

static void replace_window_text (windata_t *vwin, const char *pbuf)
{
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gtk_text_buffer_set_text(tbuf, pbuf, -1);
#else
    gtk_text_freeze(GTK_TEXT(vwin->w));
    gtk_editable_delete_text(GTK_EDITABLE(vwin->w), 0, -1);
    gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
		    NULL, NULL, pbuf, strlen(pbuf));
    gtk_text_thaw(GTK_TEXT(vwin->w));
#endif
}

static void series_view_print (windata_t *vwin)
{
    const char *pbuf;
    PRN *prn;
    series_view *sview = (series_view *) vwin->data;
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

    pbuf = gretl_print_get_buffer(prn);
    replace_window_text(vwin, pbuf);

    gretl_print_destroy(prn);
}

static void multi_series_view_print (windata_t *vwin)
{
    const char *pbuf;
    PRN *prn;
    multi_series_view *mview = (multi_series_view *) vwin->data;
    int i, vi, t, s;
    int err = 0;

    if (bufopen(&prn)) return;

    pprintf(prn, "\n%9s", " ");
    for (i=1; i<=mview->list[0]; i++) {
	vi = mview->list[i];
	if (vi >= datainfo->v) {
	    err = 1;
	    break;
	}
	pprintf(prn, "%12s", datainfo->varname[vi]);
    }

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    pputs(prn, "\n\n"); 

    for (t=0; t<mview->npoints; t++) {
	s = mview->points[t].obsnum;
	if (s >= datainfo->n) {
	    err = 1;
	    break;
	}
	print_obs_marker(s, datainfo, prn);
	for (i=1; i<=mview->list[0]; i++) {
	    vi = mview->list[i];
	    pprintf(prn, "%#12.4g", Z[vi][s]);
	}
	pputc(prn, '\n');    
    }

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    pbuf = gretl_print_get_buffer(prn);
    replace_window_text(vwin, pbuf);

    gretl_print_destroy(prn);
}

static int compare_points (const void *a, const void *b)
{
    const data_point *pa = (const data_point *) a;
    const data_point *pb = (const data_point *) b;
     
    return (pa->val > pb->val) - (pa->val < pb->val);
}

static int compare_mpoints (const void *a, const void *b)
{
    const multi_point *pa = (const multi_point *) a;
    const multi_point *pb = (const multi_point *) b;
     
    return (pa->val > pb->val) - (pa->val < pb->val);
}

static void series_view_sort (windata_t *vwin, guint action, GtkWidget *w)
{
    series_view *sview = (series_view *) vwin->data;
    
    if (series_view_allocate(sview)) {
	return;
    }

    /* sort the data */
    qsort((void *) sview->points, (size_t) sview->npoints, 
	  sizeof sview->points[0], compare_points);

    /* print sorted data to buffer */
    series_view_print(vwin);
}

void series_view_sort_by (GtkWidget *w, windata_t *vwin)
{
    multi_series_view *mview = (multi_series_view *) vwin->data;
    int v;

    if (mview == NULL || mview->list == NULL) {
	return;
    }

    if (multi_series_view_allocate(mview)) {
	return;
    }

    v = select_var_from_list(mview->list, _("Variable to sort by"));
    if (v < 0) {
	return;
    }

    mview->sortvar = v;
    mview_fill_points(mview);

    qsort((void *) mview->points, (size_t) mview->npoints, 
	  sizeof mview->points[0], compare_mpoints);

    multi_series_view_print(vwin);
}

static void 
series_view_graph (windata_t *vwin, guint action, GtkWidget *w)
{
    series_view *sview = (series_view *) vwin->data;

    if (dataset_is_time_series(datainfo)) {
	do_graph_var(sview->varnum);
    } else {
	do_boxplot_var(sview->varnum);
    }
}

static void 
scalar_to_clipboard (windata_t *vwin, guint action, GtkWidget *w)
{
    series_view *sview = (series_view *) vwin->data;
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
    series_view *sview;

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

multi_series_view *multi_series_view_new (int *list)
{
    multi_series_view *mview;

    mview = malloc(sizeof *mview);

    if (mview != NULL) {
	mview->list = list;
	mview->sortvar = 0;
	mview->npoints = 0;
	mview->points = NULL;
    } else {
	free(list);
    }

    return mview;
}

static 
void series_view_format_cancel (GtkWidget *w, series_view *sview)
{
    sview->digits = -1;
}

static 
void series_view_get_figures (GtkWidget *w, series_view *sview)
{
    sview->digits = gtk_spin_button_get_value_as_int
	(GTK_SPIN_BUTTON(sview->digit_spin));
}

static void 
set_series_float_format (GtkWidget *w, gpointer p)
{
    gint i;
    series_view *sview = (series_view *) p;

    if (GTK_TOGGLE_BUTTON(w)->active) {
        i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
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
    series_view *sview = (series_view *) vwin->data;

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

    /* spinner for number of digits */
    tmp = gtk_label_new(_("Show"));
    adj = gtk_adjustment_new(sview->digits, 1, 10, 1, 1, 1);
    sview->digit_spin = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    g_signal_connect (adj, "value_changed",
		      G_CALLBACK (series_view_get_figures), sview);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), sview->digit_spin, FALSE, FALSE, 5);

    /* select decimal places versus significant figures */
    tmp = gtk_radio_button_new_with_label (NULL, _("significant figures"));
    gtk_box_pack_start (GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    if (sview->format == 'G')
        gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (tmp), TRUE);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", GINT_TO_POINTER('G'));

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (tmp));
    tmp = gtk_radio_button_new_with_label(group, _("decimal places"));
    gtk_box_pack_start (GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    if (sview->format == 'f')
        gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (tmp), TRUE);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", GINT_TO_POINTER('f')); 

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
    gtk_window_set_transient_for(GTK_WINDOW(w), GTK_WINDOW(vwin->dialog));
#endif

    gretl_set_window_modal(w);

    gtk_main(); /* block */

    if (sview->digits > 0) {
	series_view_print(vwin);
    } else { 
	/* canceled */
	sview->digits = 6;
    }
}

