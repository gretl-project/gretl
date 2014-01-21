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

/* menustate.c: status of menus etc. */

#include "gretl.h"
#include "console.h"
#include "guiprint.h"
#include "ssheet.h"
#include "selector.h"
#include "varinfo.h"
#include "treeutils.h"
#include "session.h"
#include "gretl_ipc.h"
#include "menustate.h"

void refresh_data (void)
{
    if (data_status) {
	populate_varlist();
	set_sample_label(dataset);
    }
}

void flip (GtkUIManager *ui, const char *path, gboolean s)
{
    if (ui != NULL) {
	GtkAction *a = gtk_ui_manager_get_action(ui, path);

	if (a != NULL) {
	    gtk_action_set_sensitive(a, s);
	} else {
	    fprintf(stderr, I_("Failed to flip state of \"%s\"\n"), path);
	}
    }
}

/* by using gretl_set_window_modal() we make the main
   window visibly insensitive */

static int modcount;

static void increment_modal_count (GtkWidget *w)
{
    if (modcount == 0) {
	gtk_widget_set_sensitive(mdata->main, FALSE);
    }

    modcount++;
}

static void decrement_modal_count (GtkWidget *w, gpointer p)
{
    if (modcount > 0) {
	modcount--;
    }

    if (modcount == 0) {
	gtk_widget_set_sensitive(mdata->main, TRUE);
    }
}

void gretl_set_window_modal (GtkWidget *w)
{
    gtk_window_set_modal(GTK_WINDOW(w), TRUE);
    increment_modal_count(w);
    g_signal_connect(G_OBJECT(w), "destroy", 
		     G_CALLBACK(decrement_modal_count),
		     NULL);
}

void gretl_set_window_quasi_modal (GtkWidget *w)
{
    increment_modal_count(w);
    g_signal_connect(G_OBJECT(w), "destroy", 
		     G_CALLBACK(decrement_modal_count),
		     NULL);
}

void variable_menu_state (gboolean s)
{
    if (mdata == NULL || mdata->ui == NULL) return;

    flip(mdata->ui, "/menubar/Variable", s);
    flip(mdata->ui, "/menubar/View/xcorrgm",  
	 dataset_is_time_series(dataset));
}

static void view_items_state (gboolean s)
{
    const char *viewpaths[] = {
	"GraphVars",
	"MultiPlots",
	"summary",
	"corr",
	"xtab",
	"pca",
	"mahal",
	NULL
    };
    char fullpath[32];
    int i;

    for (i=0; viewpaths[i] != NULL; i++) {
	sprintf(fullpath, "/menubar/View/%s", viewpaths[i]);
	flip(mdata->ui, fullpath, s);
    }

    flip(mdata->ui, "/menubar/View/IconView", have_session_objects());

    flip(mdata->ui, "/menubar/View/xcorrgm",  
	 data_status && dataset_is_time_series(dataset));
}

void dataset_menubar_state (gboolean s)
{
    if (mdata == NULL || mdata->ui == NULL) return;

    flip(mdata->ui, "/menubar/File/AppendData", s);
    flip(mdata->ui, "/menubar/File/ClearData", s);
    flip(mdata->ui, "/menubar/File/SaveData", s);
    flip(mdata->ui, "/menubar/File/SaveDataAs", s);
    flip(mdata->ui, "/menubar/File/ExportData", s);
    flip(mdata->ui, "/menubar/File/MailData", s);
    flip(mdata->ui, "/menubar/Data", s);
    flip(mdata->ui, "/menubar/Add", s);
    flip(mdata->ui, "/menubar/Sample", s);
    flip(mdata->ui, "/menubar/Variable", s);
    flip(mdata->ui, "/menubar/Model", s);

    view_items_state(s);

    if (s || !have_session_objects()) {
	/* Either we're enabling dataset items, in which
	   case we should also enable the /View menu, or
	   we're disabling the dataset items and there are 
	   no session objects, in which case /View should 
	   be disabled.
	*/
	flip(mdata->ui, "/menubar/View", s);
    }

    flip(mdata->ui, "/menubar/File/NewData", !s);

    set_main_colheads_clickable(s);
}

void iconview_menubar_state (gboolean s)
{
    if (s) {
	GtkAction *a = gtk_ui_manager_get_action(mdata->ui, "/menubar/View");

	if (a != NULL && !gtk_action_get_sensitive(a)) {
	    gtk_action_set_sensitive(a, TRUE);
	}
    }

    flip(mdata->ui, "/menubar/View/IconView", s);
}

#define COMPACTABLE(d) (d->structure == TIME_SERIES && \
                        (d->pd == 4 || d->pd == 12 || \
                         d->pd == 5 || d->pd == 6 || \
                         d->pd == 7 || d->pd == 24))

#define EXPANSIBLE(d) (d->structure == TIME_SERIES && (d->pd == 1 || d->pd == 4))

#define extended_ts(d) ((d)->structure == TIME_SERIES || \
			(d)->structure == SPECIAL_TIME_SERIES || \
			(d)->structure == STACKED_TIME_SERIES)

void time_series_menu_state (gboolean s)
{
    gboolean sx = extended_ts(dataset);
    gboolean panel = dataset_is_panel(dataset);
    gboolean realpan = multi_unit_panel_sample(dataset);
    gboolean ur;

    if (mdata->ui == NULL) {
	return;
    }

    /* FIXME: we (may) need to enable/disable function
       packages that have menu attachments here.
    */

    /* unit-root tests: require time-series or panel data,
       and a time series length greater than 5
    */
    if (panel) {
	ur = dataset->pd > 5;
    } else {
	ur = s && sample_size(dataset) > 5;
    }

    /* Plots */
    flip(mdata->ui, "/menubar/View/GraphVars/TSPlot", sx);
    flip(mdata->ui, "/menubar/View/MultiPlots/MultiTS", sx);
    flip(mdata->ui, "/menubar/Variable/VarTSPlot", sx && !realpan);
    flip(mdata->ui, "/menubar/Variable/PanPlot", realpan);

    /* Variable menu */
    flip(mdata->ui, "/menubar/Variable/URTests", ur); 
    if (ur && !s) {
	/* time-series only "ur" option */
	flip(mdata->ui, "/menubar/Variable/URTests/fractint", s);
    }
    flip(mdata->ui, "/menubar/Variable/URTests/levinlin", ur && panel);
    flip(mdata->ui, "/menubar/Variable/corrgm", s);
    flip(mdata->ui, "/menubar/Variable/pergm", s);
    flip(mdata->ui, "/menubar/Variable/Filter", s);
#ifdef HAVE_X12A
    flip(mdata->ui, "/menubar/Variable/X12A", get_x12a_ok());
#endif
#ifdef HAVE_TRAMO
    flip(mdata->ui, "/menubar/Variable/Tramo", get_tramo_ok());
#endif
    flip(mdata->ui, "/menubar/Variable/Hurst", s);

    /* Model menu */
    flip(mdata->ui, "/menubar/Model/TSModels", s);

    /* Sample menu */
    flip(mdata->ui, "/menubar/Data/DataCompact", 
	 s && (COMPACTABLE(dataset) || dated_weekly_data(dataset)));
    flip(mdata->ui, "/menubar/Data/DataExpand", s && EXPANSIBLE(dataset));
}

void panel_menu_state (gboolean s)
{
    if (mdata->ui != NULL) {
	flip(mdata->ui, "/menubar/Add/AddUnit", s);
	flip(mdata->ui, "/menubar/Add/UnitDums", s);
	flip(mdata->ui, "/menubar/Add/TimeDums", s);
	flip(mdata->ui, "/menubar/Add/RangeDum", !s);
	flip(mdata->ui, "/menubar/Model/PanelModels", s);
	flip(mdata->ui, "/menubar/Model/LimdepModels/probit/reprobit", s);
	if (s && dataset->pd <= 2) {
	    flip(mdata->ui, "/menubar/Model/PanelModels/dpanel", 0);
	}
    }
}

void ts_or_panel_menu_state (gboolean s)
{
    if (mdata->ui == NULL) return;

    flip(mdata->ui, "/menubar/Data/DataSort", !s);

    flip(mdata->ui, "/menubar/Add/AddTime", s);
    flip(mdata->ui, "/menubar/Add/lags", s);
    flip(mdata->ui, "/menubar/Add/diff", s);
    flip(mdata->ui, "/menubar/Add/ldiff", s);

    s = dataset_is_seasonal(dataset);
    flip(mdata->ui, "/menubar/Add/sdiff", s);
    flip(mdata->ui, "/menubar/Add/PeriodDums", s);
}

void session_menu_state (gboolean s)
{
    if (mdata->ui != NULL) {
	flip(mdata->ui, "/menubar/View/IconView", s);
	if (!s || session_is_modified()) {
	    flip(mdata->ui, "/menubar/File/SessionFiles/SaveSession", s);
	}
	if (!s || session_is_open()) {
	    flip(mdata->ui, "/menubar/File/SessionFiles/SaveSessionAs", s);
	}
    }

    if (!s && mdata->main != NULL) {
	set_main_window_title(NULL, FALSE);
    }
}

void restore_sample_state (gboolean s)
{
    if (mdata->ui != NULL) {
	flip(mdata->ui, "/menubar/Sample/FullRange", s);
    }
}

void drop_obs_state (gboolean s)
{
    if (mdata->ui != NULL) {
	flip(mdata->ui, "/menubar/Data/RemoveObs", s);
    }
}

void compact_data_state (gboolean s)
{
    if (mdata->ui != NULL) {
	flip(mdata->ui, "/menubar/Data/DataCompact", s);
    }
}

void main_menus_enable (gboolean s)
{
    if (mdata->ui != NULL) {
	flip(mdata->ui, "/menubar/File", s);
	flip(mdata->ui, "/menubar/Tools", s);
	flip(mdata->ui, "/menubar/Data", s);
	flip(mdata->ui, "/menubar/View", s);
	flip(mdata->ui, "/menubar/Add", s);
	flip(mdata->ui, "/menubar/Sample", s);
	flip(mdata->ui, "/menubar/Variable", s);
	flip(mdata->ui, "/menubar/Model", s);
    }
}

void check_var_labels_state (GtkMenuItem *item, gpointer p)
{
    gboolean s = FALSE;

    if (dataset != NULL && dataset->v > 0) {
	if (dataset->v > 2 || strcmp(dataset->varname[1], "index")) {
	    s = TRUE;
	}
    }

    flip(mdata->ui, "/menubar/Data/VarLabels", s);
}

static gint var_popup_click (GtkWidget *w, gpointer p)
{
    gchar *item = (gchar *) p;
    int v = mdata_active_var();

    if (!strcmp(item, _("Display values"))) { 
	display_var();
    } else if (!strcmp(item, _("Summary statistics"))) { 
	do_menu_op(VAR_SUMMARY, NULL, OPT_NONE);
    } else if (!strcmp(item, _("Time series plot"))) { 
	do_graph_var(v);
    } else if (!strcmp(item, _("Panel plot..."))) { 
	do_graph_var(v);
    } else if (!strcmp(item, _("Frequency distribution"))) { 
	do_freq_dist();
    } else if (!strcmp(item, _("Boxplot"))) {
	menu_boxplot_callback(v);
    } else if (!strcmp(item, _("Gini coefficient"))) {
	do_gini();
    } else if (!strcmp(item, _("Correlogram"))) {
	do_corrgm();
    } else if (!strcmp(item, _("Periodogram"))) { 
	do_pergm(NULL);
    } else if (!strcmp(item, _("ARIMA model"))) {
	selector_set_varnum(v);
	modelspec_dialog(ARMA);
    } else if (!strcmp(item, _("Dickey-Fuller test"))) { 
	unit_root_test(ADF);
    } else if (!strcmp(item, _("KPSS test"))) { 
	unit_root_test(KPSS);
    } else if (!strcmp(item, _("Edit attributes"))) {  
	varinfo_dialog(v);
    } else if (!strcmp(item, _("Edit values"))) {  
	show_spreadsheet(SHEET_EDIT_VARLIST);
    } else if (!strcmp(item, _("Copy to clipboard"))) { 
	csv_selected_to_clipboard();
    } else if (!strcmp(item, _("Delete"))) { 
	delete_single_var(v);
    } else if (!strcmp(item, _("Add log"))) {
	add_logs_etc(LOGS, v);
    } else if (!strcmp(item, _("Add difference"))) {
	add_logs_etc(DIFF, v);
    } else if (!strcmp(item, _("Define new variable..."))) { 
	genr_callback();
    }

    gtk_widget_destroy(mdata->popup);

    return FALSE;
}

static int missvals_in_selection (void)
{
    int *list = main_window_selection_as_list();
    int miss = 0;

    if (list != NULL) {
	int i, vi, t;

	for (i=1; i<=list[0] && !miss; i++) {
	    vi = list[i];
	    for (t=dataset->t1; t<=dataset->t2; t++) {
		if (na(dataset->Z[vi][t])) {
		    miss = 1;
		    break;
		}
	    }
	}
	free(list);
    }

    return miss;
}

static int uniform_corr_option (const gchar *title, gretlopt *popt)
{
    const char *opts[] = {
	N_("Ensure uniform sample size"),
	NULL
    };
    int uniform = 0;
    int resp;

    resp = checks_only_dialog(title, NULL, opts, 1, 
			      &uniform, CORR, NULL);

    if (!canceled(resp) && uniform) {
	*popt = OPT_U;
    }

    return resp;
}

static gint selection_popup_click (GtkWidget *w, gpointer p)
{
    gchar *item = (gchar *) p;
    int ci = 0;

    if (!strcmp(item, _("Summary statistics"))) {
	ci = SUMMARY;
    } else if (!strcmp(item, _("Correlation matrix"))) {
	ci = CORR;
    }

    if (ci == CORR && missvals_in_selection()) {
	gchar *title = g_strdup_printf("gretl: %s", _("correlation matrix"));
	gretlopt opt = OPT_NONE;
	int resp;

	resp = uniform_corr_option(title, &opt);
	if (!canceled(resp)) {
	    char *buf = main_window_selection_as_string();

	    if (buf != NULL) {
		do_menu_op(ci, buf, opt);
		free(buf);
	    }
	}	    
	g_free(title);
    } else if (ci != 0) {
	char *buf = main_window_selection_as_string();
	
	if (buf != NULL) {
	    do_menu_op(ci, buf, OPT_NONE);
	    free(buf);
	}
    } else if (!strcmp(item, _("Display values"))) { 
	display_selected(); 
    } else if (!strcmp(item, _("Cross-correlogram")))  {
	xcorrgm_callback();
    } else if (!strcmp(item, _("Time series plot"))) { 
	plot_from_selection(GR_PLOT);
    } else if (!strcmp(item, _("XY scatterplot")))  {
	plot_from_selection(GR_XY);
    } else if (!strcmp(item, _("Copy to clipboard"))) { 
	csv_selected_to_clipboard();
    } else if (!strcmp(item, _("Edit values")))  {
 	show_spreadsheet(SHEET_EDIT_VARLIST);
    } else if (!strcmp(item, _("Delete")))  {
	delete_selected_vars();
    } else if (!strcmp(item, _("Add logs")))  {
	add_logs_etc(LOGS, 0);
    } else if (!strcmp(item, _("Add differences")))  {
	add_logs_etc(DIFF, 0);
    } else if (!strcmp(item, _("Define list"))) { 
	make_list_from_main();
    } else if (!strcmp(item, _("Define new variable..."))) { 
	genr_callback();
    }

    gtk_widget_destroy(mdata->popup);

    return FALSE;
}

GtkWidget *build_var_popup (void)
{
    const char *items[] = {
	N_("Display values"),
	N_("Summary statistics"),
	N_("Time series plot"),
	N_("Panel plot..."),
	N_("Frequency distribution"),
	N_("Boxplot"),
	N_("Correlogram"),
	N_("Periodogram"),
	N_("Edit attributes"),
	N_("Edit values"),
	N_("Copy to clipboard"),
	N_("Delete"),
	NULL,
	N_("Add log"),
	N_("Add difference"),
	N_("Define new variable...")
    };
    GtkWidget *menu;
    GtkWidget *item;
    int i, n = G_N_ELEMENTS(items);
    int real_panel = multi_unit_panel_sample(dataset);

    menu = gtk_menu_new();

    for (i=0; i<n; i++) {
	if (items[i] == NULL) {
	    item = gtk_separator_menu_item_new();
	    gtk_widget_show(item);
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	    continue;
	}
	if (real_panel && (i == 2 || i == 5)) {
	    /* don't offer regular ts or boxplot */
	    continue;
	}
	if (!real_panel && i == 3) {
	    /* don't offer panel plot */
	    continue;
	}	
	if ((i == 6 || i == 7) && !dataset_is_time_series(dataset)) {
	    continue;
	}
	if ((i == 2 || i == 14) && !extended_ts(dataset)) {
	    continue;
	}
	if (i == 5 && dataset_is_time_series(dataset)) {
	    /* skip boxplot option */
	    continue;
	}
	item = gtk_menu_item_new_with_label(_(items[i]));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(var_popup_click),
			 _(items[i]));
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    }

    return menu;
}

GtkWidget *build_selection_popup (void)
{
    const char *items[] = {
	N_("Display values"),
	N_("Summary statistics"),
	N_("Correlation matrix"),
	N_("Cross-correlogram"),
	N_("Time series plot"),
	N_("XY scatterplot"),
	N_("Edit values"),
	N_("Copy to clipboard"),
	N_("Delete"),
	NULL,
	N_("Add logs"),
	N_("Add differences"),
	NULL,
	N_("Define list"),
	N_("Define new variable...")
    };
    GtkWidget *menu;
    GtkWidget *item;
    int i, n = G_N_ELEMENTS(items);

    menu = gtk_menu_new();

    for (i=0; i<n; i++) {
	if (items[i] == NULL) {
	    item = gtk_separator_menu_item_new();
	    gtk_widget_show(item);
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	    continue;
	}
	if (!dataset_is_time_series(dataset) && (i == 3 || i == 11)) {
	    continue;
	}
	if (!extended_ts(dataset) && i == 4) {
	    continue;
	}
	item = gtk_menu_item_new_with_label(_(items[i]));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(selection_popup_click),
			 _(items[i]));
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    }

    return menu;
}

void clear_sample_label (void)
{
    GtkWidget *dlabel = g_object_get_data(G_OBJECT(mdata->main), "dlabel");

    gtk_label_set_text(GTK_LABEL(mdata->status), "");
    gtk_label_set_text(GTK_LABEL(dlabel), _(" No datafile loaded "));
}

void set_main_window_title (const char *name, gboolean modified)
{
#ifdef GRETL_PID_FILE
    int seqno = gretl_sequence_number();
#else
    int seqno = 0;
#endif
    gchar *title = NULL;

    if (seqno <= 1 && name == NULL) {
	gtk_window_set_title(GTK_WINDOW(mdata->main), "gretl");
    } else if (name == NULL) {
	title = g_strdup_printf("gretl (%d)", seqno);
    } else {
	gchar *prog;

	if (seqno > 1) {
	    prog = g_strdup_printf("gretl (%d)", seqno);
	} else {
	    prog = g_strdup("gretl");
	}
	
	if (!g_utf8_validate(name, -1, NULL)) {
	    gchar *trname = my_filename_to_utf8(name);
	    
	    if (modified) {
		title = g_strdup_printf("%s: %s *", prog, trname);
	    } else {
		title = g_strdup_printf("%s: %s", prog, trname);
	    }
	    g_free(trname);
	} else {
	    if (modified) {
		title = g_strdup_printf("%s: %s *", prog, name);
	    } else {
		title = g_strdup_printf("%s: %s", prog, name);
	    }
	}

	g_free(prog);
    }

    if (title != NULL) {
	gtk_window_set_title(GTK_WINDOW(mdata->main), title);
	g_free(title);
    }
}

static const char *get_pd_string (DATASET *dset)
{
    char *pdstr;

    if (custom_time_series(dset)) {
	pdstr = N_("Time series");
    } else if (dataset_is_time_series(dset)) {
	switch (dset->pd) {
	case 1:
	    pdstr = N_("Annual"); break;
	case 4:
	    pdstr = N_("Quarterly"); break;
	case 12:
	    pdstr = N_("Monthly"); break;
	case 24:
	    pdstr = N_("Hourly"); break;
	case 52:
	    pdstr = N_("Weekly"); break;
	case 5:
	case 6:
	case 7:
	    pdstr = N_("Daily"); break;
	case 10:
	    pdstr = N_("Decennial"); break;
	default:
	    pdstr = N_("Unknown"); break;
	}
    } else if (dataset_is_panel(dset)) {
	pdstr = N_("Panel");
    } else {
	pdstr = N_("Undated");
    }
    
    return pdstr;
}

void set_sample_label (DATASET *dset)
{
    GtkWidget *dlabel;
    char tmp[256];
    int tsubset;

    if (mdata == NULL) {
	return;
    }

    /* set the sensitivity of various menu items */

    time_series_menu_state(dataset_is_time_series(dset));
    panel_menu_state(dataset_is_panel(dset));
    ts_or_panel_menu_state(dataset_is_time_series(dset) ||
			   dataset_is_panel(dset));
    flip(mdata->ui, "/menubar/Data/DataTranspose", !dataset_is_panel(dset));

    tsubset = dset->t1 > 0 || dset->t2 < dset->n - 1;

    /* construct label showing summary of dataset/sample info
       (this goes at the foot of the window) 
    */

    if (complex_subsampled() && !tsubset && dataset_is_cross_section(dset)) {
	sprintf(tmp, _("Undated: Full range n = %d; current sample"
		       " n = %d"), get_full_length_n(), dataset->n);
	gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
    } else if (complex_subsampled() && dataset_is_panel(dset)) {
	char t1str[OBSLEN], t2str[OBSLEN];
	const char *pdstr = get_pd_string(dset);

	ntodate(t1str, dset->t1, dset);
	ntodate(t2str, dset->t2, dset);
	sprintf(tmp, _("%s; sample %s - %s"), _(pdstr), t1str, t2str);
	gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
    } else {
	char t1str[OBSLEN], t2str[OBSLEN];
	const char *pdstr = get_pd_string(dset);

	if (tsubset && calendar_data(dset)) {
	    /* it's too verbose to print both full range and sample */
	    ntodate(t1str, dset->t1, dset);
	    ntodate(t2str, dset->t2, dset);
	    sprintf(tmp, _("%s; sample %s - %s"), _(pdstr), t1str, t2str);
	    gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
	} else {
	    ntodate(t1str, 0, dset);
	    ntodate(t2str, dset->n - 1, dset);
	    sprintf(tmp, _("%s: Full range %s - %s"), _(pdstr), 
		    t1str, t2str);
	    if (tsubset) {
		gchar *fulltext;

		ntodate(t1str, dset->t1, dset);
		ntodate(t2str, dset->t2, dset);
		fulltext = g_strdup_printf(_("%s; sample %s - %s"), tmp, t1str, t2str);
		gtk_label_set_text(GTK_LABEL(mdata->status), fulltext);
		g_free(fulltext);
	    } else {
		gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
	    }
	}
    }

    /* construct label with datafile name (this goes above the
       data series window) */

    dlabel = g_object_get_data(G_OBJECT(mdata->main), "dlabel");

    if (dlabel != NULL) {
	if (strlen(datafile) > 2) {
	    /* data file open already */
	    const char *p = strrchr(datafile, SLASH);
	    gchar *trfname;

	    if (p != NULL) {
		trfname = my_filename_to_utf8(p + 1);
	    } else {
		trfname = my_filename_to_utf8(datafile);
	    }

	    strcpy(tmp, " ");

	    if (data_status & SESSION_DATA) {
		sprintf(tmp + 1, "Imported %s", trfname);
	    } else if (data_status & MODIFIED_DATA) {
		sprintf(tmp + 1, "%s *", trfname);
	    } else {
		sprintf(tmp + 1, "%s", trfname);
	    }

	    gtk_label_set_text(GTK_LABEL(dlabel), tmp);
	    g_free(trfname);
	} else if (data_status & MODIFIED_DATA) {
	    strcpy(tmp, _(" Unsaved data "));
	    gtk_label_set_text(GTK_LABEL(dlabel), tmp);
	}
    }

    if (complex_subsampled() || dset->t1 > 0 || 
	dset->t2 < dset->n - 1) {
	restore_sample_state(TRUE);
    }

    console_record_sample(dataset);
}

void action_entry_init (GtkActionEntry *entry)
{
    entry->stock_id = NULL;
    entry->accelerator = NULL;
    entry->tooltip = NULL;
    entry->callback = NULL;
}

int vwin_add_ui (windata_t *vwin, GtkActionEntry *entries,
		 gint n_entries, const gchar *ui_info)
{
    GtkActionGroup *actions;
    GError *err = NULL;

    actions = gtk_action_group_new("MyActions");
    gtk_action_group_set_translation_domain(actions, "gretl");

    gtk_action_group_add_actions(actions, entries, n_entries, vwin);

    vwin->ui = gtk_ui_manager_new();
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    gtk_window_add_accel_group(GTK_WINDOW(vwin->main), 
			       gtk_ui_manager_get_accel_group(vwin->ui));

    gtk_ui_manager_add_ui_from_string(vwin->ui, ui_info, -1, &err);
    if (err != NULL) {
	g_message("building menus failed: %s", err->message);
	g_error_free(err);
    }

    vwin->mbar = gtk_ui_manager_get_widget(vwin->ui, "/menubar");

    return 0;
}

static GtkActionGroup *get_named_group (GtkUIManager *uim,
					const char *name, 
					int *newgroup)
{
    GList *list = gtk_ui_manager_get_action_groups(uim); 
    GtkActionGroup *actions = NULL;

    while (list != NULL) {
	GtkActionGroup *group = list->data;

	if (!strcmp(gtk_action_group_get_name(group), name)) {
	    actions = group;
	    break;
	}
	list = list->next;
    } 

    if (actions == NULL) {
	actions = gtk_action_group_new(name);
	gtk_action_group_set_translation_domain(actions, "gretl");
	*newgroup = 1;
    } else {
	*newgroup = 0;
    }

    return actions;
}

int vwin_menu_add_item_unique (windata_t *vwin, 
			       const gchar *aname, 
			       const gchar *path, 
			       GtkActionEntry *entry)
{
    GList *list = gtk_ui_manager_get_action_groups(vwin->ui);
    GtkActionGroup *actions;
    guint id;

    while (list != NULL) {
	GtkActionGroup *group = list->data;

	if (!strcmp(aname, gtk_action_group_get_name(group))) {
	    gtk_ui_manager_remove_action_group(vwin->ui, group);
	    break;
	}
	list = list->next;
    }

    id = gtk_ui_manager_new_merge_id(vwin->ui);
    actions = gtk_action_group_new(aname);
    gtk_action_group_set_translation_domain(actions, "gretl");

    gtk_action_group_add_actions(actions, entry, 1, vwin);
    gtk_ui_manager_add_ui(vwin->ui, id, path, entry->name, entry->name,
			  GTK_UI_MANAGER_MENUITEM, FALSE);

    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);
	
    return id;
}

/* Retrieve existing "AdHoc" action group from @uim, or add a 
   new group of this name to the UIManager and return it.
*/

GtkActionGroup *get_ad_hoc_group (GtkUIManager *uim,
				  int *newgroup)
{
    return get_named_group(uim, "AdHoc", newgroup);
}

/* Adds the specified @entry to vwin->ui at @path; returns
   the "merge_id", which can be used to remove the item
*/

int vwin_menu_add_item (windata_t *vwin, const gchar *path, 
			GtkActionEntry *entry)
{
    GtkActionGroup *actions;
    int newgroup = 1;
    guint id;

    actions = get_ad_hoc_group(vwin->ui, &newgroup);
    gtk_action_group_add_actions(actions, entry, 1, vwin);
    id = gtk_ui_manager_new_merge_id(vwin->ui); 

    gtk_ui_manager_add_ui(vwin->ui, id, path, entry->name, entry->name,
			  GTK_UI_MANAGER_MENUITEM, FALSE);

    if (newgroup) {
	gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
	g_object_unref(actions);
    }

    return id;
}

/* Adds the specified @entries to vwin->ui at @path; returns
   the "merge_id", which can be used to remove the items
*/

int vwin_menu_add_items (windata_t *vwin, const gchar *path, 
			 GtkActionEntry *entries, int n)
{
    GtkActionGroup *actions;
    int newgroup = 1;
    guint id;
    int i;

    actions = get_ad_hoc_group(vwin->ui, &newgroup);
    gtk_action_group_add_actions(actions, entries, n, vwin);
    id = gtk_ui_manager_new_merge_id(vwin->ui);

    for (i=0; i<n; i++) {
	gtk_ui_manager_add_ui(vwin->ui, id, path, 
			      entries[i].name, entries[i].name,
			      GTK_UI_MANAGER_MENUITEM, FALSE);
    }

    if (newgroup) {
	gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
	g_object_unref(actions);
    }

    return id;
}

int vwin_menu_add_radios (windata_t *vwin, const gchar *path, 
			  GtkRadioActionEntry *entries, int n,
			  int deflt, GCallback callback)
{
    guint id = gtk_ui_manager_new_merge_id(vwin->ui);
    GtkActionGroup *actions;
    int i;

    actions = gtk_action_group_new("Radios");
    gtk_action_group_set_translation_domain(actions, "gretl");

    gtk_action_group_add_radio_actions(actions, entries, n,
				       deflt, callback,
				       vwin);
    for (i=0; i<n; i++) {
	gtk_ui_manager_add_ui(vwin->ui, id, path, 
			      entries[i].name, entries[i].name,
			      GTK_UI_MANAGER_MENUITEM, FALSE);
    }
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    return id;
}

int vwin_menu_add_menu (windata_t *vwin, const gchar *path, 
			GtkActionEntry *entry)
{
    guint id = gtk_ui_manager_new_merge_id(vwin->ui);
    GtkActionGroup *actions;
    gchar *grpname;
    static int seq;

    grpname = g_strdup_printf("NewMenu%d", seq);
    actions = gtk_action_group_new(grpname);
    gtk_action_group_set_translation_domain(actions, "gretl");
    g_free(grpname);
    seq++;

    gtk_action_group_add_actions(actions, entry, 1, vwin);
    gtk_ui_manager_add_ui(vwin->ui, id, path, entry->name, entry->name,
			  GTK_UI_MANAGER_MENU, FALSE);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    return id;
}

void vwin_menu_add_separator (windata_t *vwin, const gchar *path)
{
    guint id = gtk_ui_manager_new_merge_id(vwin->ui);

    gtk_ui_manager_add_ui(vwin->ui, id, path, NULL, NULL,
			  GTK_UI_MANAGER_SEPARATOR, FALSE);
}
