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

/* series_view.c for gretl */

#include "gretl.h"
#include "textutil.h"
#include "dlgutils.h"
#include "menustate.h"
#include "textbuf.h"
#include "series_view.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#else
# include "clipboard.h"
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
    int err = 0;

    if (sview->npoints != 0) {
	/* already allocated */
	return 0;
    } else {
	int t, tp, T = datainfo->t2 - datainfo->t1 + 1;
	int v = sview->varnum;

	sview->points = malloc(T * sizeof *sview->points);
	if (sview->points == NULL) {
	    err = E_ALLOC;
	} else {
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
    }

    return err;
}

static int multi_series_view_allocate (multi_series_view *mview)
{
    if (mview->npoints != 0) {
	/* already allocated */
	return 0;
    } else {
	int T = datainfo->t2 - datainfo->t1 + 1;

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

static PRN *series_view_print_csv (windata_t *vwin)
{
    series_view *sview = (series_view *) vwin->data;
    char dchar = datainfo->delim;
    PRN *prn;
    int t;

    if (bufopen(&prn)) {
	return NULL;
    }
    
    pprintf(prn, "obs%c%s\n", dchar, datainfo->varname[sview->varnum]);
    for (t=0; t<sview->npoints; t++) {
	if (na(sview->points[t].val)) {
	    pprintf(prn, "%s%cNA\n", sview->points[t].label, dchar);
	} else {
	    pprintf(prn, "\"%s\"%c%.8g\n", sview->points[t].label, dchar,
		    sview->points[t].val);
	}
    }

    return prn;
}

/* for printing sorted or reformatted data, window with a
   single series */

static void series_view_print (windata_t *vwin)
{
    series_view *sview = (series_view *) vwin->data;
    const char *pbuf;
    int obslen = 8;
    PRN *prn;
    int t, len;

    if (bufopen(&prn)) return;

    for (t=0; t<sview->npoints; t++) {
	len = strlen(sview->points[t].label);
	if (len > obslen) {
	    obslen = len;
	}
    }
	
    /* print formatted data to buffer */

    if (1) { /* FIXME scalars */
	pprintf(prn, "\n%*s ", obslen, _("Obs"));
	pprintf(prn, "%13s\n\n", datainfo->varname[sview->varnum]);
	for (t=0; t<sview->npoints; t++) {
	    if (na(sview->points[t].val)) {
		pprintf(prn, "%*s\n", obslen, sview->points[t].label);
	    } else if (sview->format == 'G') {
		pprintf(prn, "%*s %#13.*g\n", obslen,
			sview->points[t].label,
			sview->digits, sview->points[t].val);
	    } else {
		pprintf(prn, "%*s %13.*f\n", obslen,
			sview->points[t].label,
			sview->digits, sview->points[t].val);
	    }
	}
    } else {
	if (sview->format == 'G') {
	    pprintf(prn, "\n%s = %#13.*g", 
		    datainfo->varname[sview->varnum], 
		    sview->digits, Z[sview->varnum][0]);
	} else {
	    pprintf(prn, "\n%s = %13.*f", 
		    datainfo->varname[sview->varnum], 
		    sview->digits, Z[sview->varnum][0]);
	}
    }

    pbuf = gretl_print_get_buffer(prn);
    textview_set_text(vwin->text, pbuf);

    gretl_print_destroy(prn);
}

static int *make_obsvec (multi_series_view *mview)
{
    int *ov;
    int t;

    ov = mymalloc((mview->npoints + 1) * sizeof *ov);
    if (ov != NULL) {
	ov[0] = mview->npoints;
	for (t=0; t<mview->npoints; t++) {
	    ov[t+1] = mview->points[t].obsnum;
	}
    }

    return ov;
}

static void multi_series_view_print_sorted (windata_t *vwin)
{
    multi_series_view *mview = (multi_series_view *) vwin->data;
    const char *pbuf;
    int *obsvec;
    PRN *prn;
    int err = 0;

    obsvec = make_obsvec(mview);
    if (obsvec == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	free(obsvec);
	return;
    }

    err = print_data_sorted(mview->list, obsvec, (const double **) Z, 
			    datainfo, prn);
    if (err) {
	gui_errmsg(err);
    } else {
	pbuf = gretl_print_get_buffer(prn);
	textview_set_text(vwin->text, pbuf);
    }

    free(obsvec);
    gretl_print_destroy(prn);
}

int series_view_is_sorted (windata_t *vwin)
{
    multi_series_view *mview = (multi_series_view *) vwin->data;

    return mview->sortvar != 0;
}

PRN *vwin_print_sorted_as_csv (windata_t *vwin)
{
    multi_series_view *mview;
    int *obsvec;
    PRN *prn;
    int err = 0;

    if (vwin->role == VIEW_SERIES) {
	return series_view_print_csv(vwin);
    }

    mview = (multi_series_view *) vwin->data;

    obsvec = make_obsvec(mview);
    if (obsvec == NULL) {
	return NULL;
    }

    if (bufopen(&prn)) {
	free(obsvec);
	return NULL;
    }

    gretl_print_set_format(prn, GRETL_FORMAT_CSV);
    err = print_data_sorted(mview->list, obsvec, (const double **) Z, 
			    datainfo, prn);
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	prn = NULL;
    } 

    free(obsvec);

    return prn;
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

void series_view_sort (GtkWidget *w, windata_t *vwin)
{
    series_view *sview = (series_view *) vwin->data;
    
    if (series_view_allocate(sview)) {
	return;
    }

    qsort(sview->points, sview->npoints, 
	  sizeof sview->points[0], compare_points);

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

    qsort(mview->points, mview->npoints, 
	  sizeof mview->points[0], compare_mpoints);

    multi_series_view_print_sorted(vwin);
}

void series_view_graph (GtkWidget *w, windata_t *vwin)
{
    series_view *sview = (series_view *) vwin->data;

    if (dataset_is_time_series(datainfo)) {
	do_graph_var(sview->varnum);
    } else {
	do_boxplot_var(sview->varnum);
    }
}

void scalar_to_clipboard (windata_t *vwin)
{
    series_view *sview = (series_view *) vwin->data;
    double val;
    gchar *buf;

    val = Z[sview->varnum][0];

    if (sview->format == 'G') {
	buf = g_strdup_printf("%#.*g", sview->digits, val);
    } else {
	buf = g_strdup_printf("%.*f", sview->digits, val);
    }

#ifdef G_OS_WIN32
    win_buf_to_clipboard(buf);
#else
    buf_to_clipboard(buf);
#endif

    g_free(buf);
}

int *series_view_get_list (windata_t *vwin)
{
    int *list = NULL;

    if (vwin->role == VIEW_SERIES) {
	series_view *sview = vwin->data;

	if (sview != NULL) {
	    list = gretl_list_new(1);
	    if (list != NULL) {
		list[1] = sview->varnum;
	    }
	}
    } else {
	multi_series_view *mview = (multi_series_view *) vwin->data;

	if (mview != NULL) {
	    list = gretl_list_copy(mview->list);
	}
    }

    return list;
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

int has_sortable_data (windata_t *vwin)
{
    multi_series_view *mview;

    if (vwin == NULL || vwin->role != PRINT || vwin->data == NULL) {
	return 0;
    }

    mview = vwin->data;

    return mview->list != NULL && mview->list[0] <= 5;
}

multi_series_view *multi_series_view_new (const int *list)
{
    multi_series_view *mview = NULL;

    if (list == NULL) {
	return NULL;
    } 

    mview = malloc(sizeof *mview);

    if (mview != NULL) {
	mview->list = gretl_list_copy(list);
	if (mview->list == NULL) {
	    free(mview);
	    mview = NULL;
	} else {
	    mview->sortvar = 0;
	    mview->npoints = 0;
	    mview->points = NULL;
	}
    } 

    return mview;
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

void series_view_format_dialog (GtkWidget *src, windata_t *vwin)
{
    GtkWidget *dlg;
    GtkWidget *tmp, *hbox; 
    GtkObject *adj;
    GSList *group;
    series_view *sview = (series_view *) vwin->data;

    if (series_view_allocate(sview)) {
	return;
    }

    dlg = gretl_dialog_new(_("gretl: data format"), NULL,
			   GRETL_DLG_BLOCK | GRETL_DLG_MODAL);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Select data format"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show_all(hbox);

    /* spinner for number of digits */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Show"));
    adj = gtk_adjustment_new(sview->digits, 1, 10, 1, 1, 1);
    sview->digit_spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    g_signal_connect (adj, "value_changed",
		      G_CALLBACK(series_view_get_figures), sview);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), sview->digit_spin, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show_all(hbox);

    /* select decimal places versus significant figures */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_radio_button_new_with_label(NULL, _("significant figures"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    if (sview->format == 'G') {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
    }
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", GINT_TO_POINTER('G'));
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, TRUE, TRUE, 0);
    gtk_widget_show_all(hbox);

    hbox = gtk_hbox_new(FALSE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    tmp = gtk_radio_button_new_with_label(group, _("decimal places"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    if (sview->format == 'f') {
        gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (tmp), TRUE);
    }
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_series_float_format), sview);
    g_object_set_data(G_OBJECT(tmp), "action", GINT_TO_POINTER('f'));
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, TRUE, TRUE, 0);
    gtk_widget_show_all(hbox);

    /* Cancel button */
    cancel_options_button(GTK_DIALOG(dlg)->action_area, dlg, 
			  &sview->digits);
   
    /* OK button */
    tmp = ok_button(GTK_DIALOG(dlg)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    gtk_widget_show(dlg);

    if (sview->digits > 0) {
	series_view_print(vwin);
    } else { 
	/* canceled */
	sview->digits = 6;
    }
}

