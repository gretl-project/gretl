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
    data_point_t *points;
} series_view_t;

/* .................................................................. */

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

void build_series_view_popup (windata_t *win)
{
    if (win->popup != NULL) return;

    win->popup = gtk_menu_new();

    add_popup_item(_("Sort values"), win->popup, 
		   G_CALLBACK(series_view_sort), 
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
	vwin->data = sview;
    }
}
