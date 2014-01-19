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

/* progress bar implementation for gretl */

#include <gtk/gtk.h>

#ifdef UPDATER
# include <stdio.h>
# include <stdlib.h>
# include "updater.h"
#else
# include "libgretl.h"
#endif

typedef struct _ProgressData {
    GtkWidget *window;
    GtkWidget *label;
    GtkWidget *pbar;
    int *cancel;
} ProgressData;

static void destroy_progress (GtkWidget *widget, ProgressData **ppdata)
{
    (*ppdata)->window = NULL;
    free(*ppdata);
    *ppdata = NULL;
}

#ifdef UPDATER

static void cancel_progress (GtkWidget *widget, ProgressData *pdata)
{
    *pdata->cancel = 1;
    gtk_widget_destroy(pdata->window);
}

#endif

static ProgressData *build_progress_window (int flag, int *cancel)
{
    ProgressData *pdata;
    GtkWidget *align, *vbox;
#ifdef UPDATER
    GtkWidget *separator;
    GtkWidget *button;
#endif

    pdata = malloc(sizeof *pdata);
    if (pdata == NULL) {
	return NULL;
    }

    pdata->cancel = cancel;

    pdata->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_resizable(GTK_WINDOW(pdata->window), FALSE);

    if (flag == SP_LOAD_INIT) {
	gtk_window_set_title(GTK_WINDOW(pdata->window), _("gretl: loading data"));
    } else if (flag == SP_SAVE_INIT) {
	gtk_window_set_title(GTK_WINDOW(pdata->window), _("gretl: storing data"));
    } else if (flag == SP_FONT_INIT) {
	gtk_window_set_title(GTK_WINDOW(pdata->window), _("gretl: scanning fonts"));
    }
	
    gtk_container_set_border_width(GTK_CONTAINER(pdata->window), 0);

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_container_add(GTK_CONTAINER(pdata->window), vbox);

    /* Add a label */
    pdata->label = gtk_label_new("");
    gtk_box_pack_start(GTK_BOX(vbox), pdata->label, FALSE, FALSE, 0);
        
    /* Create a centering alignment object */
    align = gtk_alignment_new(0.5, 0.5, 0, 0);
    gtk_box_pack_start(GTK_BOX(vbox), align, FALSE, FALSE, 5);

    /* Create the GtkProgressBar */
    pdata->pbar = gtk_progress_bar_new();
    gtk_container_add(GTK_CONTAINER(align), pdata->pbar);

    /* Add separator and cancel button? */
#ifdef UPDATER
    separator = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), separator, FALSE, FALSE, 0);

    button = gtk_button_new_with_label(_("Cancel"));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(cancel_progress),
		     pdata);
    gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);

    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_widget_grab_default(button);
#endif

    gtk_widget_show_all(pdata->window);

    return pdata;
}

int show_progress (gint64 res, gint64 expected, int flag)
{
    static ProgressData *pdata;
    static gint64 offs;
    static int cancel;

    if (expected == 0) {
	return SP_RETURN_DONE;
    }

    if (res < 0 || flag == SP_FINISH) {
	/* clean up and get out */
	if (pdata != NULL && pdata->window != NULL) {
	    gtk_widget_destroy(GTK_WIDGET(pdata->window)); 
	    while (gtk_events_pending()) {
		gtk_main_iteration();
	    }	    
	}
	return SP_RETURN_DONE;
    }

    if (flag == SP_LOAD_INIT || flag == SP_SAVE_INIT || flag == SP_FONT_INIT) {
	/* initialize the progress bar */
	gchar *bytestr = NULL;

	offs = 0;
	cancel = 0;

	pdata = build_progress_window(flag, &cancel);
	if (pdata == NULL) {
	    return 0; 
	}

	g_signal_connect(G_OBJECT(pdata->window), "destroy",
			 G_CALLBACK(destroy_progress),
			 &pdata);

	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(pdata->pbar), 0.0);

	if (flag == SP_LOAD_INIT) {
#ifdef G_OS_WIN32
	    bytestr = g_strdup_printf("%s %I64d Kbytes", _("Retrieving"),
				      expected / 1024);
#else
	    bytestr = g_strdup_printf("%s %lld Kbytes", _("Retrieving"),
				      (long long) expected / 1024);
#endif
	} else if (flag == SP_SAVE_INIT) {
#ifdef G_OS_WIN32
	    bytestr = g_strdup_printf("%s %I64d Kbytes", _("Storing"),
				      expected / 1024);
#else
	    bytestr = g_strdup_printf("%s %lld Kbytes", _("Storing"),
				      (long long) expected / 1024);
#endif
	} else if (flag == SP_FONT_INIT) {
	    bytestr = g_strdup_printf(_("Scanning %d fonts"), (int) expected);
	}

	gtk_label_set_text(GTK_LABEL(pdata->label), bytestr);
	g_free(bytestr);

	while (gtk_events_pending()) {
	    gtk_main_iteration();
	}
    }

    if ((flag == SP_NONE || flag == SP_TOTAL) && cancel) {
	/* the user canceled */
	cancel = 0;
	return SP_RETURN_CANCELED;
    }    

    if ((flag == SP_NONE || flag == SP_TOTAL) &&
	(pdata == NULL || pdata->window == NULL)) {
	/* something has gone wrong */
	return 0;
    }

    if (flag == SP_TOTAL) {
	offs = res;
    } else {
	offs += res;
    }

    if (offs > expected && pdata != NULL) {
	gtk_widget_destroy(GTK_WIDGET(pdata->window)); 
	return SP_RETURN_DONE;
    }

    if (offs <= expected && pdata != NULL) {
	/* display the completed fraction */
	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(pdata->pbar), 
				      (gdouble) ((double) offs / expected));
	while (gtk_events_pending()) {
	    gtk_main_iteration();
	}
    } else {
	if (pdata != NULL && pdata->window != NULL) {
	    gtk_widget_destroy(GTK_WIDGET(pdata->window)); 
	}
	return SP_RETURN_DONE;
    }
	
    return SP_RETURN_OK;
}

