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
#include "uservar.h"
#include "treeutils.h"
#include "session.h"
#include "gretl_ipc.h"
#include "gretl_www.h"
#include "fncall.h"
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
            fprintf(stderr, "Failed to flip state of \"%s\"\n", path);
        }
    }
}

void menu_item_set_tooltip (GtkUIManager *ui, const char *path,
			    const char *tip)
{
    if (ui != NULL) {
	GtkWidget *w;

	w = gtk_ui_manager_get_widget(ui, path);
	if (w != NULL) {
	    gtk_widget_set_tooltip_text(w, _(tip));
	}
    }
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
        "Summary",
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

static void gfn_menuitems_state (void)
{
    GtkActionGroup *ag;
    GList *aglist;
    DataReq dreq;
    int err;

    aglist = gtk_ui_manager_get_action_groups(mdata->ui);

    while (aglist != NULL) {
        ag = aglist->data;
        if (GPOINTER_TO_INT(g_object_get_data(G_OBJECT(ag), "datachk"))) {
            dreq = pkg_get_data_requirement(ag);
            err = check_function_needs(dataset, dreq, 0, NULL);
            gtk_action_group_set_sensitive(ag, !err);
	    if (err) {
		gretl_error_clear();
	    }
        }
        aglist = aglist->next;
    }
}

void dataset_menubar_state (gboolean s)
{
    static int mail_ok = -1;

    if (mdata == NULL || mdata->ui == NULL) return;

    if (mail_ok < 0) {
        mail_ok = curl_does_smtp();
    }

    flip(mdata->ui, "/menubar/File/AppendData", s);
    flip(mdata->ui, "/menubar/File/ClearData", s);
    flip(mdata->ui, "/menubar/File/SaveData", s);
    flip(mdata->ui, "/menubar/File/SaveDataAs", s);
    flip(mdata->ui, "/menubar/File/ExportData", s);
    flip(mdata->ui, "/menubar/File/MailData", mail_ok && s);
    flip(mdata->ui, "/menubar/Data", s);
    flip(mdata->ui, "/menubar/Add", s);
    flip(mdata->ui, "/menubar/Sample", s);
    flip(mdata->ui, "/menubar/Variable", s);
    flip(mdata->ui, "/menubar/Model", s);

    view_items_state(s);
    gfn_menuitems_state();

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

    if (s) {
	time_series_menu_state(dataset_is_time_series(dataset));
    }

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

#define OK_MIDAS_PD(p) (p == 1 || p == 4 || p == 12)

#define COMPACTABLE(d) ((d->structure == TIME_SERIES || \
                         d->structure == SPECIAL_TIME_SERIES) && \
                        d->pd > 1 && d->pd != 10)

#define EXPANSIBLE(d) (d->structure == TIME_SERIES && \
                       (d->pd == 1 || d->pd == 4))

#define extended_ts(d) ((d)->structure == TIME_SERIES || \
                        (d)->structure == SPECIAL_TIME_SERIES || \
                        (d)->structure == STACKED_TIME_SERIES)

void time_series_menu_state (gboolean s)
{
    gboolean sx = extended_ts(dataset);
    gboolean panel = dataset_is_panel(dataset);
    gboolean realpan = multi_unit_panel_sample(dataset);
    gboolean have_map = dataset_get_mapfile(dataset) != NULL;
    gboolean ur;

    if (mdata->ui == NULL) {
        return;
    }

    /* enable/disable function packages that have menu
       attachments */
    gfn_menuitems_state();

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
    flip(mdata->ui, "/menubar/View/GraphVars/MapPlot", have_map);

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
    flip(mdata->ui, "/menubar/Variable/BDS", s);
    flip(mdata->ui, "/menubar/Variable/tdisagg", s &&
         quarterly_or_monthly(dataset));

    /* Model menu */
    flip(mdata->ui, "/menubar/Model/TSModels", s);
    flip(mdata->ui, "/menubar/Model/TSMulti", s);
    flip(mdata->ui, "/menubar/Model/TSModels/midasreg",
         s && OK_MIDAS_PD(dataset->pd));

    /* Data menu */
    flip(mdata->ui, "/menubar/Data/DataCompact",
         s && (COMPACTABLE(dataset) || dated_weekly_data(dataset)));
    flip(mdata->ui, "/menubar/Data/DataExpand", s && EXPANSIBLE(dataset));
}

void panel_menu_state (gboolean s)
{
    if (mdata->ui != NULL) {
        flip(mdata->ui, "/menubar/Add/Panel", s);
        flip(mdata->ui, "/menubar/Add/RangeDum", !s);
        flip(mdata->ui, "/menubar/Add/Time/idxvals", !s);
        flip(mdata->ui, "/menubar/Add/Time/PeriodDums", !s);
        flip(mdata->ui, "/menubar/Model/PanelModels", s);
        flip(mdata->ui, "/menubar/Model/LimdepModels/probit/reprobit", s);
        flip(mdata->ui, "/menubar/Model/PanelModels/dpanel",
             s && dataset->pd > 2);
        gfn_menuitems_state();
    }
}

void ts_or_panel_menu_state (gboolean s)
{
    if (mdata->ui == NULL) return;

    flip(mdata->ui, "/menubar/Data/DataSort", !s);
    flip(mdata->ui, "/menubar/Add/Time", s);
    flip(mdata->ui, "/menubar/Add/Time/idxvals",
         s && !dataset_is_panel(dataset));

    s = dataset_is_seasonal(dataset);
    if (!s && dataset_is_seasonal_panel(dataset)) {
        s = dataset->panel_pd == 4 ||
            dataset->panel_pd == 12 ||
            dataset->panel_pd == 24;
    }

    flip(mdata->ui, "/menubar/Add/Time/sdiff", s);
    flip(mdata->ui, "/menubar/Add/Time/PeriodDums", s);
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

void sample_related_menu_state (void)
{
    if (mdata->ui != NULL) {
	gboolean s = dataset_is_subsampled(dataset);

        /* actions supported only when subsampled */
	flip(mdata->ui, "/menubar/Sample/FullRange", s);
	flip(mdata->ui, "/menubar/Sample/PermaSample", s);
        /* actions supported when not subsampled */
	flip(mdata->ui, "/menubar/Data/DataStructure", !s);
        flip(mdata->ui, "/menubar/Data/DataTranspose", !s);

	if (s && dataset->submask != NULL && dataset->submask == RESAMPLED) {
	    flip(mdata->ui, "/menubar/Sample/PermaSample", FALSE);
	}
	if (!s && dataset_is_panel(dataset)) {
	    flip(mdata->ui, "/menubar/Data/DataTranspose", FALSE);
	}
    }
}

void sample_menubar_state (gboolean s)
{
    if (mdata->ui != NULL) {
        flip(mdata->ui, "/menubar/Sample", s);
    }
}

void drop_obs_state (gboolean s)
{
    if (mdata->ui != NULL) {
        flip(mdata->ui, "/menubar/Data/RemoveObs", s);
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

static int missvals_in_selection (const int *list)
{
    int i, vi, t;
    int ret = 0;

    for (i=1; i<=list[0] && !ret; i++) {
        vi = list[i];
        for (t=dataset->t1; t<=dataset->t2; t++) {
            if (na(dataset->Z[vi][t])) {
                ret = 1;
                break;
            }
        }
    }

    return ret;
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
        *popt = OPT_N;
    }

    return resp;
}

static void right_click_corr (void)
{
    int *list = main_window_selection_as_list();
    gretlopt opt = OPT_NONE;
    char *buf;

    if (list != NULL && list[0] > 2 && missvals_in_selection(list)) {
        gchar *title;
        int resp;

        title = g_strdup_printf("gretl: %s", _("correlation matrix"));
        resp = uniform_corr_option(title, &opt);
        g_free(title);

        if (canceled(resp)) {
            return;
        }
    }

    free(list);
    buf = main_window_selection_as_string();

    if (buf != NULL) {
        do_menu_op(CORR, buf, opt, NULL);
        free(buf);
    }
}

static int is_integer_valued (int v)
{
    const double *x = dataset->Z[v];
    double x0 = NADBL;
    int nonconst = 0;
    int t, n = 0;

    for (t=dataset->t1; t<=dataset->t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (!ok_int(x[t])) {
            /* out of integer bounds */
            return 0;
        }
        if (x[t] != floor(x[t])) {
            /* non-integral */
            return 0;
        }
        if (na(x0)) {
            x0 = x[t];
        } else if (x[t] != x0) {
            nonconst = 1;
        }
        n++;
    }

    return n > 2 && nonconst;
}

int series_is_dummifiable (int v)
{
    if (series_is_discrete(dataset, v)) {
        /* must be OK */
        return 1;
    } else if (gretl_isdummy(0, dataset->n-1, dataset->Z[v])) {
        /* already a 0/1 dummy */
        return 0;
    } else {
        /* could be OK? */
        return is_integer_valued(v);
    }
}

static int dataset_could_be_midas (const DATASET *dset)
{
    if (dataset_is_time_series(dset) &&
        (dset->pd == 1 || dset->pd == 4 || dset->pd == 12)) {
        return 1;
    } else {
        return 0;
    }
}

/* regular context menu items */

enum MenuIdx_ {
    MNU_MPLOT,
    MNU_MSAVE,
    MNU_DISP,
    MNU_EDIT,
    MNU_STATS,
    MNU_TPLOT,
    MNU_PPLOT,
    MNU_FDIST,
    MNU_BPLOT,
    MNU_CGRAM,
    MNU_PGRAM,
    MNU_ATTRS,
    MNU_STRS,
    MNU_CORR,
    MNU_COND,
    MNU_XCORR,
    MNU_SCATR,
    MNU_CLIPB,
    MNU_DELET,
    MNU_SEPAR,
    MNU_LOGS,
    MNU_DIFF,
    MNU_PCDIF,
    MNU_TMEAN,
    MNU_XMEAN,
    MNU_IDXV,
    MNU_DUMIF,
    MNU_GENR,
    MNU_LIST,
    MNU_TDIS
};

/* MIDAS-special context menu items */

enum MDSIdx_ {
    MDS_DISP,
    MDS_TPLOT,
    MDS_LOGS,
    MDS_DIFF,
    MDS_SEPAR,
    MDS_CDISP,
    MDS_CPLOT,
    MDS_CEDIT,
    MDS_CDEL,
    MDS_GENR,
    MDS_LIST
};

enum MenuTarg_ {
    T_SINGLE, /* item applicable only for a single series */
    T_MULTI,  /* item applicable only for multiple series */
    T_BOTH,   /* item applicable in both cases */
};

typedef enum MenuIdx_ MenuIdx;
typedef enum MDSIdx_ MDSIdx;
typedef enum MenuTarg_ MenuTarg;

struct popup_entries {
    MenuIdx idx;       /* one of the MenuIdxvalues above */
    const char *str;   /* translatable string */
    MenuTarg target;   /* one of the MenuTarget values above */
};

struct mpopup_entries {
    MDSIdx idx;        /* one of the MDSIdx values above */
    const char *str;   /* translatable string */
};

struct popup_entries main_pop_entries[] = {
    { MNU_MPLOT, N_("Display map..."), T_BOTH },
    { MNU_MSAVE, N_("Write map as GeoJSON..."), T_BOTH },
    { MNU_SEPAR, NULL, T_BOTH },
    { MNU_DISP,  N_("Display values"), T_BOTH },
    { MNU_EDIT,  N_("Edit values"), T_BOTH },
    { MNU_STATS, N_("Summary statistics"), T_BOTH },
    { MNU_TPLOT, N_("Time series plot"), T_BOTH },
    { MNU_PPLOT, N_("Panel plot..."), T_SINGLE },
    { MNU_FDIST, N_("Frequency distribution"), T_SINGLE },
    { MNU_BPLOT, N_("Boxplot"), T_SINGLE },
    { MNU_CGRAM, N_("Correlogram"), T_SINGLE, },
    { MNU_PGRAM, N_("Periodogram"), T_SINGLE, },
    { MNU_ATTRS, N_("Edit attributes"), T_SINGLE },
    { MNU_STRS,  N_("Show string table"), T_SINGLE },
    { MNU_CORR,  N_("Correlation matrix"), T_MULTI },
    { MNU_COND,  N_("Collinearity"), T_MULTI },
    { MNU_XCORR, N_("Cross-correlogram"), T_MULTI },
    { MNU_SCATR, N_("XY scatterplot"), T_MULTI },
    { MNU_CLIPB, N_("Copy to clipboard"), T_BOTH },
    { MNU_DELET, N_("Delete"), T_BOTH },
    { MNU_SEPAR, NULL, T_BOTH },
    { MNU_LOGS,  N_("Add log"), T_SINGLE },
    { MNU_DIFF,  N_("Add difference"), T_SINGLE },
    { MNU_PCDIF, N_("Add percent change..."), T_SINGLE },
    { MNU_TMEAN, N_("Add time mean..."), T_SINGLE },
    { MNU_XMEAN, N_("Add cross-sectional mean..."), T_SINGLE },
    { MNU_IDXV,  N_("Add index values..."), T_SINGLE },
    { MNU_DUMIF, N_("Dummify..."), T_SINGLE },
    { MNU_LOGS,  N_("Add logs"), T_MULTI },
    { MNU_DIFF,  N_("Add differences"), T_MULTI },
    { MNU_PCDIF, N_("Add percent changes..."), T_MULTI },
    { MNU_IDXV,  N_("Add index values..."), T_MULTI },
    { MNU_TDIS,  N_("Disaggregate..."), T_SINGLE },
    { MNU_SEPAR, NULL, T_BOTH },
    { MNU_GENR,  N_("Define new variable..."), T_BOTH },
    { MNU_LIST,  N_("Define list"), T_MULTI }
};

struct mpopup_entries midas_pop_entries[] = {
    { MDS_DISP,  N_("Display values") },
    { MDS_TPLOT, N_("Time series plot") },
    { MDS_LOGS,  N_("Add logs...") },
    { MDS_DIFF,  N_("Add differences...") },
    { MDS_SEPAR, NULL },
    { MDS_CDISP, N_("Display components") },
    { MDS_CPLOT, N_("Plot components") },
    { MDS_CEDIT, N_("Edit components") },
    { MDS_CDEL,  N_("Delete components") },
    { MDS_SEPAR, NULL },
    { MDS_GENR,  N_("Define new variable...") },
    { MDS_LIST,  N_("Define list") }
};

static gint var_popup_click (GtkWidget *w, gpointer p)
{
    MenuIdx i = GPOINTER_TO_INT(p);
    int v = mdata_active_var();

    switch (i) {
    case MNU_DISP:
        display_var();
        break;
    case MNU_STATS:
        do_menu_op(VAR_SUMMARY, NULL, OPT_NONE, NULL);
        break;
    case MNU_TPLOT:
    case MNU_PPLOT:
        do_graph_var(v);
        break;
    case MNU_MPLOT:
        geoplot_callback();
        break;
    case MNU_MSAVE:
        map_save_callback();
        break;
    case MNU_FDIST:
        do_freq_dist();
        break;
    case MNU_BPLOT:
        menu_boxplot_callback(v);
        break;
    case MNU_CGRAM:
        do_corrgm();
        break;
    case MNU_PGRAM:
        do_pergm(NULL);
        break;
    case MNU_ATTRS:
        varinfo_dialog(v);
        break;
    case MNU_STRS:
	display_string_table(v);
	break;
    case MNU_EDIT:
        show_spreadsheet(SHEET_EDIT_VARLIST);
        break;
    case MNU_CLIPB:
        selected_series_to_clipboard();
        break;
    case MNU_DELET:
        delete_single_var(v);
        break;
    case MNU_LOGS:
    case MNU_DIFF:
        add_logs_etc(i == MNU_LOGS ? LOGS : DIFF, v, 0);
        break;
    case MNU_PCDIF:
        single_percent_change_dialog(v, 0);
        break;
    case MNU_TMEAN:
        panel_mean_dialog(v, 0);
        break;
    case MNU_XMEAN:
        panel_mean_dialog(v, 1);
        break;
    case MNU_IDXV:
        single_percent_change_dialog(v, 1);
        break;
    case MNU_DUMIF:
        add_discrete_dummies(v);
        break;
    case MNU_TDIS:
	tdisagg_dialog(v);
	break;
    case MNU_GENR:
        genr_callback();
        break;
    default:
        break;
    }

    gtk_widget_destroy(mdata->popup);

    return FALSE;
}

GtkWidget *build_var_popup (int selvar)
{
    GtkWidget *menu, *item;
    int i, j, n = G_N_ELEMENTS(main_pop_entries);
    int real_panel = multi_unit_panel_sample(dataset);
    int have_map = dataset_get_mapfile(dataset) != NULL;
    int strvar = is_string_valued(dataset, selvar);
    int nullbak = 0;

    menu = gtk_menu_new();

    for (j=0; j<n; j++) {
        if (main_pop_entries[j].target == T_MULTI) {
            /* not applicable */
            continue;
        }
        i = main_pop_entries[j].idx;
        if (i == MNU_SEPAR) {
            if (!nullbak) {
                /* don't insert two consecutive separators */
                item = gtk_separator_menu_item_new();
                gtk_widget_show(item);
                gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
                nullbak = 1;
            }
            continue;
        }
        if (real_panel && (i == MNU_TPLOT || i == MNU_BPLOT)) {
            /* don't offer regular ts or boxplot */
            continue;
        }
        if (!real_panel && (i == MNU_PPLOT ||
                            i == MNU_TMEAN ||
                            i == MNU_XMEAN)) {
            /* don't offer panel plot or panel means */
            continue;
        }
        if (!have_map && (i == MNU_MPLOT || i == MNU_MSAVE)) {
            /* don't offer map plot or save */
            continue;
        }
        if ((i == MNU_CGRAM || i == MNU_PGRAM || i == MNU_IDXV) &&
            !dataset_is_time_series(dataset)) {
            /* correlogram, periodogram, index values */
            continue;
        }
        if ((i == MNU_TPLOT || i == MNU_DIFF || i == MNU_PCDIF) &&
            !extended_ts(dataset)) {
            /* time-series plot, difference, percent change */
            continue;
        }
        if (i == MNU_BPLOT && dataset_is_time_series(dataset)) {
            /* skip boxplot option */
            continue;
        }
        if (i == MNU_DUMIF && !series_is_dummifiable(selvar)) {
            /* skip dummify option */
            continue;
        }
        if (i == MNU_STATS && strvar) {
            /* skip (numerical) summary stats option */
            continue;
        }
        if (i == MNU_STRS && !strvar) {
            /* no string table */
            continue;
        }
        if (i == MNU_TDIS && series_get_orig_pd(dataset, selvar) == 0) {
            /* skip temporal disaggregation option */
            continue;
        }
        item = gtk_menu_item_new_with_label(_(main_pop_entries[j].str));
        g_signal_connect(G_OBJECT(item), "activate",
                         G_CALLBACK(var_popup_click),
                         GINT_TO_POINTER(i));
        gtk_widget_show(item);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);

        nullbak = 0;
    }

    return menu;
}

static gint selection_popup_click (GtkWidget *w, gpointer p)
{
    MenuIdx i = GPOINTER_TO_INT(p);
    int ci = 0;

    if (i == MNU_STATS) {
        ci = POP_SUMMARY;
    } else if (i == MNU_CORR) {
        ci = CORR;
    }

    if (ci == CORR) {
        right_click_corr();
    } else if (ci != 0) {
        char *buf = main_window_selection_as_string();

        if (buf != NULL) {
            do_menu_op(ci, buf, OPT_NONE, NULL);
            free(buf);
        }
    } else if (i == MNU_DISP) {
        display_selected();
    } else if (i == MNU_COND) {
        cond_number_callback();
    } else if (i == MNU_XCORR)  {
        xcorrgm_callback();
    } else if (i == MNU_TPLOT) {
        plot_from_selection(GR_PLOT);
    } else if (i == MNU_SCATR)  {
        plot_from_selection(GR_XY);
    } else if (i == MNU_MPLOT) {
        geoplot_callback();
    } else if (i == MNU_MSAVE) {
	map_save_callback();
    } else if (i == MNU_EDIT)  {
        show_spreadsheet(SHEET_EDIT_VARLIST);
    } else if (i == MNU_CLIPB) {
        selected_series_to_clipboard();
    } else if (i == MNU_DELET)  {
        delete_selected_vars();
    } else if (i == MNU_LOGS || i == MNU_DIFF)  {
        add_logs_etc(i == MNU_LOGS ? LOGS : DIFF, 0, 0);
    } else if (i == MNU_PCDIF) {
        multi_percent_change_dialog(0);
    } else if (i == MNU_IDXV) {
        multi_percent_change_dialog(1);
    } else if (i == MNU_LIST) {
        make_list_from_main();
    } else if (i == MNU_GENR) {
        genr_callback();
    }

    gtk_widget_destroy(mdata->popup);

    return FALSE;
}

static GtkWidget *build_regular_selection_popup (void)
{
    GtkWidget *menu, *item;
    int i, j, n = G_N_ELEMENTS(main_pop_entries);
    int nullbak = 0;

    menu = gtk_menu_new();

    for (j=0; j<n; j++) {
        if (main_pop_entries[j].target == T_SINGLE) {
            /* for single selection only */
            continue;
        }
        i = main_pop_entries[j].idx;
        if (i == MNU_SEPAR) {
            if (!nullbak) {
                /* don't insert two consecutive separators */
                item = gtk_separator_menu_item_new();
                gtk_widget_show(item);
                gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
                nullbak = 1;
            }
            continue;
        }
        if ((i == MNU_TPLOT || i == MNU_XCORR) &&
            !dataset_is_time_series(dataset)) {
            continue;
        }
        if (i == MNU_TPLOT && !extended_ts(dataset)) {
            continue;
        }
        if ((i == MNU_MPLOT || i == MNU_MSAVE) && dataset->mapfile == NULL) {
            continue;
        }
        item = gtk_menu_item_new_with_label(_(main_pop_entries[j].str));
        g_signal_connect(G_OBJECT(item), "activate",
                         G_CALLBACK(selection_popup_click),
                         GINT_TO_POINTER(i));
        gtk_widget_show(item);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
        nullbak = 0;
    }

    return menu;
}

static gint midas_popup_click (GtkWidget *w, gpointer p)
{
    MDSIdx i = GPOINTER_TO_INT(p);

    if (i == MDS_DISP || i == MDS_TPLOT) {
        int *list = main_window_selection_as_list();

        midas_list_callback(list, NULL, i == MDS_DISP ? PRINT : PLOT);
        free(list);
    } else if (i == MDS_LOGS || i == MDS_DIFF)  {
        add_logs_etc(i == MDS_LOGS ? LOGS : DIFF, 0, 1);
    } else if (i == MDS_CDISP) {
        display_selected();
    } else if (i == MDS_CPLOT) {
        plot_from_selection(GR_PLOT);
    } else if (i == MDS_CEDIT) {
        show_spreadsheet(SHEET_EDIT_VARLIST);
    } else if (i == MDS_CDEL) {
        delete_selected_vars();
    } else if (i == MDS_LIST) {
        make_list_from_main();
    } else if (i == MDS_GENR) {
        genr_callback();
    }

    gtk_widget_destroy(mdata->popup);

    return FALSE;
}

static GtkWidget *build_midas_popup (void)
{
    GtkWidget *menu, *item;
    int n = G_N_ELEMENTS(midas_pop_entries);
    int i, j;

    menu = gtk_menu_new();

    for (j=0; j<n; j++) {
        i = midas_pop_entries[j].idx;
        if (i == MDS_SEPAR) {
            item = gtk_separator_menu_item_new();
            gtk_widget_show(item);
            gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
            continue;
        }
        item = gtk_menu_item_new_with_label(_(midas_pop_entries[j].str));
        g_signal_connect(G_OBJECT(item), "activate",
                         G_CALLBACK(midas_popup_click),
                         GINT_TO_POINTER(i));
        gtk_widget_show(item);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    }

    return menu;
}

GtkWidget *build_selection_popup (void)
{
    int midas_list = 0;

    if (dataset_could_be_midas(dataset)) {
        int *list = main_window_selection_as_list();

        if (gretl_is_midas_list(list, dataset)) {
            midas_list = 1;
        }
        free(list);
    }

    if (midas_list) {
        return build_midas_popup();
    } else {
        return build_regular_selection_popup();
    }
}

void clear_sample_label (void)
{
    GtkWidget *dlabel = g_object_get_data(G_OBJECT(mdata->main), "dlabel");

    gtk_label_set_text(GTK_LABEL(mdata->status), "");
    gtk_label_set_text(GTK_LABEL(dlabel), _(" No datafile loaded "));
}

/* Note: if @name is not NULL here it will be the name of
   a gretl session file */

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
        /* show session name */
        const char *mod = modified ? " *" : "";
        gchar *prog;

        if (seqno > 1) {
            prog = g_strdup_printf("gretl (%d)", seqno);
        } else {
            prog = g_strdup("gretl");
        }
        title = g_strdup_printf("%s: session %s%s", prog, name, mod);
        g_free(prog);
    }

    if (title != NULL) {
        gtk_window_set_title(GTK_WINDOW(mdata->main), title);
        g_free(title);
    }
}

static const char *daily_pdstr (DATASET *dset)
{
    const char *pdstrs[] = {
        N_("Daily (5 days)"),
        N_("Daily (6 days)"),
        N_("Daily (7 days)"),
        N_("Daily (5 days, incomplete)"),
        N_("Daily (6 days, incomplete)"),
        N_("Daily (7 days, incomplete)"),
        N_("Daily (5 days, undated)"),
        N_("Daily (6 days, undated)"),
        N_("Daily (7 days, undated)")
    };
    int i = dset->pd - 5;
    int j = 0;

    if (dset->markers == DAILY_DATE_STRINGS) {
        /* incomplete */
        j = 1;
    } else if (dset->sd0 < 100000) {
        /* undated */
        j = 2;
    }

    return pdstrs[i + j*3];
}

static const char *get_pd_string (DATASET *dset)
{
    const char *pdstr;

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
            pdstr = daily_pdstr(dset); break;
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

static void special_panel_label (DATASET *dset)
{
    DATASET tset = {0};
    char st1[OBSLEN], st2[OBSLEN];
    int u1, u2;
    gchar *s;

    tset.structure = TIME_SERIES;
    tset.pd = dset->panel_pd;
    tset.sd0 = dset->panel_sd0;
    u1 = 1 + dset->t1 / dset->pd;
    u2 = 1 + dset->t2 / dset->pd;
    ntolabel(st1, 0, &tset);
    ntolabel(st2, dset->pd - 1, &tset);
    s = g_strdup_printf(_("Panel: units %d:%d, time %s:%s"),
			u1, u2, st1, st2);
    gtk_label_set_text(GTK_LABEL(mdata->status), s);
    g_free(s);
}

void set_sample_label (DATASET *dset)
{
    GtkWidget *dlabel;
    gchar *tmp = NULL;
    int tsubset;

    if (mdata == NULL) {
        return;
    }

    /* set the sensitivity of various menu items */

    time_series_menu_state(dataset_is_time_series(dset));
    panel_menu_state(dataset_is_panel(dset));
    ts_or_panel_menu_state(dataset_is_time_series(dset) ||
                           dataset_is_panel(dset));

    tsubset = dset->t1 > 0 || dset->t2 < dset->n - 1;

    /* construct label showing summary of dataset/sample info
       (this goes at the foot of the window)
    */

    if (dataset_is_panel(dset) && dset->panel_pd > 0) {
	special_panel_label(dset);
    } else if (complex_subsampled() && !tsubset && dataset_is_cross_section(dset)) {
        tmp = g_strdup_printf(_("Undated: Full range n = %d; current sample"
                                " n = %d"), get_full_length_n(), dataset->n);
        gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
    } else if (complex_subsampled() && dataset_is_panel(dset)) {
        char t1str[OBSLEN], t2str[OBSLEN];
        const char *pdstr = get_pd_string(dset);

        ntolabel(t1str, dset->t1, dset);
        ntolabel(t2str, dset->t2, dset);
        tmp = g_strdup_printf(_("%s; sample %s - %s"), _(pdstr), t1str, t2str);
        gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
    } else {
        char t1str[OBSLEN], t2str[OBSLEN];
        const char *pdstr = get_pd_string(dset);

        if (calendar_data(dset) && tsubset) {
            /* it's too verbose to print both full range and sample */
            ntolabel(t1str, dset->t1, dset);
            ntolabel(t2str, dset->t2, dset);
            tmp = g_strdup_printf(_("%s; sample %s - %s"), _(pdstr), t1str, t2str);
            gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
        } else if (calendar_data(dset) && complex_subsampled()) {
            /* ditto, too verbose */
            tmp = g_strdup_printf(_("%s; sample %s - %s"), _(pdstr), dset->stobs,
                                  dset->endobs);
            gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
        } else {
            int done = 0;

            ntolabel(t1str, 0, dset);
            ntolabel(t2str, dset->n - 1, dset);
            tmp = g_strdup_printf(_("%s: Full range %s - %s"), _(pdstr),
                                  t1str, t2str);
            if (dataset_is_panel(dset) && !tsubset) {
                GString *full = g_string_new(tmp);

                g_string_append(full, " ");
                g_string_append(full, _("(unit:period)"));
                gtk_label_set_text(GTK_LABEL(mdata->status), full->str);
                g_string_free(full, TRUE);
                done = 1;
            }
            if (tsubset) {
                gchar *full;

                ntolabel(t1str, dset->t1, dset);
                ntolabel(t2str, dset->t2, dset);
                full = g_strdup_printf(_("%s; sample %s - %s"), tmp, t1str, t2str);
                gtk_label_set_text(GTK_LABEL(mdata->status), full);
                g_free(full);
            } else if (!done) {
                gtk_label_set_text(GTK_LABEL(mdata->status), tmp);
            }
        }
    }

    g_free(tmp);

    /* construct label with datafile name (this goes above the
       data series window) */

    dlabel = g_object_get_data(G_OBJECT(mdata->main), "dlabel");

    if (dlabel != NULL) {
        GString *dl = NULL;

        if (strlen(datafile) > 2) {
            /* data file open already */
            gchar *basename = g_path_get_basename(datafile);

            dl = g_string_new(" ");
            if (data_status & SESSION_DATA) {
                g_string_append_printf(dl, _("Imported %s"), basename);
            } else if (data_status & MODIFIED_DATA) {
                g_string_append_printf(dl, "%s *", basename);
            } else {
                g_string_append(dl, basename);
            }
            gtk_label_set_text(GTK_LABEL(dlabel), dl->str);
            g_free(basename);
        } else if (data_status & MODIFIED_DATA) {
            dl = g_string_new(_(" Unsaved data "));
            gtk_label_set_text(GTK_LABEL(dlabel), dl->str);
        }
	if (dl != NULL) {
	    g_string_free(dl, TRUE);
	}
    }

    sample_related_menu_state();
    console_record_sample(dataset);
}

void set_workdir_label (void)
{
    GtkWidget *wlabel;

    wlabel = g_object_get_data(G_OBJECT(mdata->main), "wlabel");

    if (wlabel != NULL) {
        gchar *fmt, *wdir, *buf;
        int maxlen = swallow ? 32 : 56;

        fmt = g_strdup_printf("<span color=\"%s\">%%s</span>",
                              blue_for_text());
        wdir = g_strdup(gretl_workdir());
        trim_slash(wdir);
        if (g_utf8_strlen(wdir, -1) > maxlen) {
            gretl_utf8_truncate(wdir, maxlen - 3);
            strcat(wdir, "...");
        }
        buf = g_markup_printf_escaped(fmt, wdir);
        gtk_label_set_markup(GTK_LABEL(wlabel), buf);
        g_free(buf);
        g_free(wdir);
        g_free(fmt);
    }
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
