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

/* menustate.c: status of menus etc. */

#include "gretl.h"
#include "console.h"
#include "guiprint.h"

void refresh_data (void)
{
    if (data_status) {
	populate_varlist();
    }
}

void edit_info_state (gboolean s)
{
    flip(mdata->ifac, "/Data/Edit info", s);
}

void add_remove_markers_state (gboolean s)
{
    flip(mdata->ifac, "/Sample/Add case markers...", !s);
    flip(mdata->ifac, "/Sample/Remove case markers", s);
}

void main_menubar_state (gboolean s)
{
    if (mdata->ifac == NULL) return;

    flip(mdata->ifac, "/File/Append data", s);
    flip(mdata->ifac, "/File/Clear data set", s);
    flip(mdata->ifac, "/File/Save data", s);
    flip(mdata->ifac, "/File/Save data as", s);
    flip(mdata->ifac, "/File/Export data", s);
    flip(mdata->ifac, "/File/Create data set", !s);
    flip(mdata->ifac, "/Data", s);
    flip(mdata->ifac, "/Sample", s);
    flip(mdata->ifac, "/Variable", s);
    flip(mdata->ifac, "/Model", s);

    if (s) {
	edit_info_state(!(data_status & BOOK_DATA));
	add_remove_markers_state(datainfo->S != NULL);
    }
}

#ifndef OLD_GTK

static GtkItemFactoryEntry time_series_model_items[] = {
    { N_("/Model/Time series/_Cochrane-Orcutt..."), NULL, model_callback, CORC, NULL, GNULL },
    { N_("/Model/Time series/_Hildreth-Lu..."), NULL, model_callback, HILU, NULL, GNULL },
    { N_("/Model/Time series/_Prais-Winsten..."), NULL, model_callback, PWE, NULL, GNULL },
    { N_("/Model/Time series/_Autoregressive estimation..."), NULL, model_callback, AR, NULL, GNULL },
    { N_("/Model/Time series/ARMA_X..."), NULL, model_callback, ARMA, NULL, GNULL },
    { N_("/Model/Time series/GARCH..."), NULL, model_callback, GARCH, NULL, GNULL },
    { N_("/Model/Time series/_Vector Autoregression..."), NULL, model_callback, VAR, NULL, GNULL },
    { N_("/Model/Time series/Cointegration test"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Model/Time series/Cointegration test/Engle-Granger..."), NULL, 
      selector_callback, COINT, NULL, GNULL },
    { N_("/Model/Time series/Cointegration test/Johansen..."), NULL, 
      selector_callback, COINT2, NULL, GNULL }
};

#else

static GtkItemFactoryEntry time_series_model_items[] = {
    { N_("/Model/Time series/_Cochrane-Orcutt..."), NULL, model_callback, CORC, NULL },
    { N_("/Model/Time series/_Hildreth-Lu..."), NULL, model_callback, HILU, NULL },
    { N_("/Model/Time series/_Prais-Winsten..."), NULL, model_callback, PWE, NULL },
    { N_("/Model/Time series/_Autoregressive estimation..."), NULL, model_callback, AR, NULL },
    { N_("/Model/Time series/ARMA_X..."), NULL, model_callback, ARMA, NULL },
    { N_("/Model/Time series/GARCH..."), NULL, model_callback, GARCH, NULL },
    { N_("/Model/Time series/_Vector Autoregression..."), NULL, model_callback, VAR, NULL },
    { N_("/Model/Time series/Cointegration test"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model/Time series/Cointegration test/Engle-Granger..."), NULL, 
      selector_callback, COINT, NULL },
    { N_("/Model/Time series/Cointegration test/Johansen..."), NULL, 
      selector_callback, COINT2, NULL }
};

#endif

#define DATASET_COMPACTABLE(d) (d->structure == TIME_SERIES && \
                                (d->pd == 4 || d->pd == 12 || \
                                 d->pd == 5 || d->pd == 7))

void time_series_menu_state (gboolean s)
{
    if (mdata->ifac == NULL) {
	return;
    }

    /* Data menu */
    flip(mdata->ifac, "/Data/Graph specified vars/Time series plot...", s);

    /* Variable menu */
    flip(mdata->ifac, "/Variable/Time series plot", s);
    flip(mdata->ifac, "/Variable/Correlogram", s);
    flip(mdata->ifac, "/Variable/Spectrum", s);
    flip(mdata->ifac, "/Variable/Runs test", s);
    flip(mdata->ifac, "/Variable/Augmented Dickey-Fuller test", s);
    flip(mdata->ifac, "/Variable/KPSS test", s);
    flip(mdata->ifac, "/Variable/ARMA model", s);
#ifdef HAVE_X12A
    flip(mdata->ifac, "/Variable/X-12-ARIMA analysis", s);
#endif
#ifdef HAVE_TRAMO
    flip(mdata->ifac, "/Variable/TRAMO analysis", s);
#endif
    flip(mdata->ifac, "/Variable/Hurst exponent", s);
    /* Model menu */
    flip(mdata->ifac, "/Model/Time series", s);
    /* Sample menu */
    flip(mdata->ifac, "/Sample/Compact data...", 
	 s && DATASET_COMPACTABLE(datainfo));

    if (s) {
	GtkWidget *w =  
	    gtk_item_factory_get_widget(mdata->ifac, 
					"/Model/Time series/Cochrane-Orcutt...");

	if (w == NULL) {
	    int i, n = sizeof time_series_model_items / 
		sizeof time_series_model_items[0];

	    for (i=0; i<n; i++) {
		gtk_item_factory_create_item(mdata->ifac, 
					     &time_series_model_items[i], 
					     mdata, 1);
	    }
	}
    }
}

void panel_menu_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Model/Pooled OLS (panel)...", s);
	flip(mdata->ifac, "/Data/Add variables/unit dummies", s);
	flip(mdata->ifac, "/Data/Add variables/panel dummies", s);
    }
}

void ts_or_panel_menu_state (gboolean s)
{
    if (mdata->ifac == NULL) return;

    flip(mdata->ifac, "/Data/Add variables/time trend", s);
    flip(mdata->ifac, "/Data/Add variables/lags of selected variables", s);
    flip(mdata->ifac, "/Data/Add variables/first differences of selected variables", s);
    flip(mdata->ifac, "/Data/Add variables/log differences of selected variables", s);
    flip(mdata->ifac, "/Data/Add variables/periodic dummies", s);
}

void session_menu_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Session/Icon view", s);
	flip(mdata->ifac, "/Session/Save", s);
	flip(mdata->ifac, "/Session/Save as...", s);
    }	
}

void restore_sample_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Sample/Restore full range", s);
    }
}

void compact_data_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Sample/Compact data...", s);
    }
}

static gint var_popup_click (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (!strcmp(item, _("Display values"))) 
	display_var();
    if (!strcmp(item, _("Descriptive statistics"))) 
	do_menu_op(NULL, VAR_SUMMARY, NULL);
    else if (!strcmp(item, _("Time series plot"))) 
	do_graph_var(mdata->active_var);
    else if (!strcmp(item, _("Frequency distribution"))) 
	do_menu_op(NULL, FREQ, NULL);
    else if (!strcmp(item, _("Frequency plot"))) 
	do_freqplot(NULL, 0, NULL);
    else if (!strcmp(item, _("Boxplot")))
	do_boxplot_var(mdata->active_var);
    else if (!strcmp(item, _("Correlogram"))) 
	gretl_callback(NULL, CORRGM, NULL);
    else if (!strcmp(item, _("Spectrum"))) 
	do_pergm(NULL, 0, NULL);
    else if (!strcmp(item, _("Spectrum (Bartlett)"))) 
	do_pergm(NULL, 1, NULL);
    else if (!strcmp(item, _("ARMA model"))) 
	arma_options_dialog(NULL, 0, NULL);
    else if (!strcmp(item, _("Dickey-Fuller test"))) 
	unit_root_test(NULL, ADF, NULL);
    else if (!strcmp(item, _("KPSS test"))) 
	unit_root_test(NULL, KPSS, NULL);
    else if (!strcmp(item, _("Runs test"))) 
	do_menu_op(NULL, RUNS, NULL);
    else if (!strcmp(item, _("Hurst exponent"))) 
	do_hurst(NULL, 0, NULL);
    else if (!strcmp(item, _("Edit attributes")))  
	varinfo_dialog(mdata->active_var, 1);
    else if (!strcmp(item, _("Delete"))) 
	delete_selected_vars(mdata->active_var);
    else if (!strcmp(item, _("Define new variable..."))) 
	gretl_callback(NULL, GENR, NULL);

    gtk_widget_destroy(mdata->popup);

    return FALSE;
}

GtkWidget *build_var_popup (void)
{
    const char *var_items[]={
	N_("Display values"),
	N_("Descriptive statistics"),
	N_("Time series plot"),
	N_("Frequency distribution"),
	N_("Frequency plot"),
	N_("Boxplot"),
	N_("Correlogram"),
	N_("Spectrum"),
	N_("Spectrum (Bartlett)"),
	N_("ARMA model"),
	N_("Dickey-Fuller test"),
	N_("KPSS test"),
	N_("Runs test"),
	N_("Hurst exponent"),
	N_("Edit attributes"),
	N_("Delete"),
	N_("Define new variable...")
    };

    GtkWidget *var_menu;
    GtkWidget *var_item;
    int i, n_items = sizeof var_items / sizeof var_items[0];

    var_menu = gtk_menu_new();

    for (i=0; i<n_items; i++) {
	if (!dataset_is_time_series(datainfo) && 
	    (i == 2 || (i >= 6 && i <= 11))) {
	    continue;
	}
	var_item = gtk_menu_item_new_with_label(_(var_items[i]));
	g_signal_connect(G_OBJECT(var_item), "activate",
			 G_CALLBACK(var_popup_click),
			 _(var_items[i]));
	gtk_widget_show(var_item);
#ifndef OLD_GTK
	gtk_menu_shell_append(GTK_MENU_SHELL(var_menu), var_item);
#else
	GTK_WIDGET_SET_FLAGS (var_item, GTK_SENSITIVE | GTK_CAN_FOCUS); /* ?? */
	gtk_menu_append(GTK_MENU(var_menu), var_item);
#endif
    }

    return var_menu;
}

static gint selection_popup_click (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (!strcmp(item, _("Display values"))) 
	display_selected(NULL, 0, NULL); 
    else if (!strcmp(item, _("Descriptive statistics"))) 
	do_menu_op(NULL, SUMMARY_SELECTED, NULL);
    else if (!strcmp(item, _("Correlation matrix"))) 
	do_menu_op(NULL, CORR_SELECTED, NULL);
    else if (!strcmp(item, _("Principal components"))) 
	do_menu_op(NULL, PCA, NULL);
    else if (!strcmp(item, _("Mahalanobis distances"))) 
	do_menu_op(NULL, MAHAL, NULL);
    else if (!strcmp(item, _("Time series plot"))) 
	plot_from_selection(NULL, GR_PLOT, NULL);
    else if (!strcmp(item, _("XY scatterplot"))) 
	plot_from_selection(NULL, GR_XY, NULL);
    else if (!strcmp(item, _("Copy to clipboard"))) 
	csv_selected_to_clipboard();
    else if (!strcmp(item, _("Delete"))) 
	delete_selected_vars(0);

    gtk_widget_destroy(mdata->popup);

    return FALSE;
}

GtkWidget *build_selection_popup (void)
{
    const char *items[] = {
	N_("Display values"),
	N_("Descriptive statistics"),
	N_("Correlation matrix"),
	N_("Time series plot"),
	N_("XY scatterplot"),
	N_("Principal components"),
	N_("Mahalanobis distances"),
	N_("Copy to clipboard"),
	N_("Delete"),
    };

    GtkWidget *sel_menu;
    GtkWidget *item;
    int i, n_items = sizeof items / sizeof items[0];

    sel_menu = gtk_menu_new();

    for (i=0; i<n_items; i++) {
	if (!dataset_is_time_series(datainfo) && i == 3) {
	    continue;
	}
	item = gtk_menu_item_new_with_label(_(items[i]));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(selection_popup_click),
			 _(items[i]));
	gtk_widget_show(item);
#ifndef OLD_GTK
	gtk_menu_shell_append(GTK_MENU_SHELL(sel_menu), item);
#else 
	gtk_menu_append(GTK_MENU(sel_menu), item);
#endif
    }

    return sel_menu;
}

/* ......................................................... */

void clear_sample_label (void)
{
    GtkWidget *dlabel = g_object_get_data(G_OBJECT(mdata->w), "dlabel");

    gtk_label_set_text(GTK_LABEL(mdata->status), "");
    gtk_label_set_text(GTK_LABEL(dlabel), _(" No datafile loaded "));
}

void set_sample_label (DATAINFO *pdinfo)
{
    GtkWidget *dlabel = g_object_get_data(G_OBJECT(mdata->w), "dlabel");
    char pdstr[16];
    char stobs[OBSLEN], endobs[OBSLEN];
    char labeltxt[80];

    ntodate(stobs, 0, pdinfo);
    ntodate(endobs, pdinfo->n - 1, pdinfo);

    if (custom_time_series(pdinfo)) {
	strcpy(pdstr, _("Time series"));
    } else if (dataset_is_time_series(pdinfo)) {
	switch (pdinfo->pd) {
	case 1:
	    strcpy(pdstr, _("Annual")); break;
	case 4:
	    strcpy(pdstr, _("Quarterly")); break;
	case 12:
	    strcpy(pdstr, _("Monthly")); break;
	case 24:
	    strcpy(pdstr, _("Hourly")); break;
	case 52:
	    strcpy(pdstr, _("Weekly")); break;
	case 5:
	case 6:
	case 7:
	    strcpy(pdstr, _("Daily")); break;
	case 10:
	    strcpy(pdstr, _("Decennial")); break;
	default:
	    strcpy(pdstr, _("Unknown")); break;
	}
    } else if (dataset_is_panel(pdinfo)) {
	strcpy(pdstr, _("Panel"));
    } else {
	strcpy(pdstr, _("Undated"));
    }

    time_series_menu_state(dataset_is_time_series(pdinfo));
    panel_menu_state(dataset_is_panel(pdinfo));
    ts_or_panel_menu_state(dataset_is_time_series(pdinfo) ||
			   dataset_is_panel(pdinfo));

    flip(mdata->ifac, "/Sample/Restructure panel...", 
	 dataset_is_panel(pdinfo));

    flip(mdata->ifac, "/Sample/Transpose data...", 
	 !dataset_is_panel(pdinfo));

    sprintf(labeltxt, _("%s: Full range %s - %s"), 
	    pdstr, stobs, endobs);

    if (pdinfo->t1 > 0 || pdinfo->t2 < pdinfo->n - 1) {
	char t1str[OBSLEN], t2str[OBSLEN];
	char biglabel[128];

	ntodate(t1str, pdinfo->t1, pdinfo);
	ntodate(t2str, pdinfo->t2, pdinfo);
	sprintf(biglabel, _("%s; sample %s - %s"), labeltxt, t1str, t2str);
	gtk_label_set_text(GTK_LABEL(mdata->status), biglabel);
    } else {
	gtk_label_set_text(GTK_LABEL(mdata->status), labeltxt);
    }

    if (strlen(paths.datfile) > 2) {
	/* data file open already */
	if (strrchr(paths.datfile, SLASH) == NULL) {
	    sprintf(labeltxt, " %s ", paths.datfile);
	} else {
	    sprintf(labeltxt, " %s ", strrchr(paths.datfile, SLASH) + 1);
	}
	if (data_status & MODIFIED_DATA) { 
	    strcat(labeltxt, "* ");
	} else if (data_status & GZIPPED_DATA) {
	    strcat(labeltxt, "[gz] ");
	}
	if (dlabel != NULL) {
	    gtk_label_set_text(GTK_LABEL(dlabel), labeltxt);
	}
    } else if (data_status & MODIFIED_DATA) {
	strcpy(labeltxt, _(" Unsaved data "));
	gtk_label_set_text(GTK_LABEL(dlabel), labeltxt);
    }

    console_record_sample(datainfo);
}




