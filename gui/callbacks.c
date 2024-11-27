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

/* callbacks.c for gretl */

#include "gretl.h"
#include "selector.h"
#include "session.h"
#include "database.h"
#include "datafiles.h"
#include "textbuf.h"
#include "textutil.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "ssheet.h"
#include "treeutils.h"
#include "datawiz.h"
#include "winstack.h"
#include "fncall.h"
#include "matrix_extra.h"

static void doubleclick_action (windata_t *vwin)
{
    switch (vwin->role) {
    case MAINWIN:
        if (dataset != NULL && dataset->n > 0) {
            display_var();
        }
        break;
    case TEXTBOOK_DATA:
        browser_open_data(NULL, vwin);
        break;
    case PS_FILES:
        browser_open_ps(NULL, vwin);
        break;
    case FUNC_FILES:
        browser_call_func(NULL, vwin);
        break;
    case NATIVE_DB:
        open_db_index(NULL, vwin);
        break;
    case REMOTE_DB:
        open_remote_db_index(NULL, vwin);
        break;
    case NATIVE_SERIES:
    case RATS_SERIES:
    case PCGIVE_SERIES:
    case REMOTE_SERIES:
        display_db_series(vwin);
        break;
    case DBNOMICS_TOP:
        open_dbnomics_provider(NULL, vwin);
        break;
    case DBNOMICS_DB:
        open_dbnomics_dataset(NULL, vwin);
        break;
    case DBNOMICS_SERIES:
        open_dbnomics_series(NULL, vwin);
        break;
    default:
        break;
    }
}

void listbox_select_row (GtkTreeSelection *selection, gpointer data)
{
    windata_t *win = (windata_t *) data;
    GtkTreeIter iter;
    GtkTreeModel *model;
    GtkTreePath *path;

    if (!gtk_tree_selection_get_selected(selection, &model, &iter)) {
        return;
    }

    path = gtk_tree_model_get_path(model, &iter);
    win->active_var = tree_path_get_row_number(path);
    gtk_tree_path_free(path);
}

gint listbox_double_click (GtkWidget *widget, GdkEventButton *event,
                           windata_t *win)
{
    if (event != NULL && event->type == GDK_2BUTTON_PRESS
        && event->button == 1) {
        doubleclick_action(win);
    }
    return FALSE;
}

gboolean listbox_drag (GtkWidget *listbox, GdkEventMotion *event,
                       gpointer data)
{
    gint x, y;
    GdkModifierType state;
    GtkTreeView *view = GTK_TREE_VIEW(listbox);
    GtkTreePath *path;

    if (event->is_hint) {
        gdk_window_get_pointer(event->window, &x, &y, &state);
    } else {
        x = event->x;
        y = event->y;
        state = event->state;
    }

    if ((state & GDK_BUTTON1_MASK) &&
        gtk_tree_view_get_path_at_pos(view, x, y, &path,
                                      NULL, NULL, NULL)) {
        GtkTreeSelection *select = NULL;
        GtkTreePath *anchor_path = NULL;
        gchar *anchor_id = NULL;
        gint row;
        int anchor;
        static gint lastrow;

        anchor = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(listbox),
                                                   "active_row"));
        row = tree_path_get_row_number(path);

        select = gtk_tree_view_get_selection(view);
        if (select == NULL) {
            return FALSE;
        }
        anchor_id = g_strdup_printf("%d", anchor);
        anchor_path = gtk_tree_path_new_from_string(anchor_id);
        g_free(anchor_id);

        if (row != lastrow) {
            gtk_tree_selection_unselect_all(select);
            gtk_tree_selection_select_range(select, anchor_path,
                                            path);
        }

        gtk_tree_path_free(path);
        gtk_tree_path_free(anchor_path);

        lastrow = row;
    }

    return FALSE;
}

struct open_data_code {
    int c;
    const gchar *s;
};

struct open_data_code open_data_codes[] = {
    { OPEN_DATA,       "OpenData" },
    { APPEND_DATA,     "AppendData" },
    { OPEN_MARKERS,    "AddMarkers" },
    { OPEN_RATS_DB,    "RATSDB" },
    { OPEN_PCGIVE_DB,  "PcGiveDB" },
    { 0, NULL }
};

static int open_data_code (const gchar *s)
{
    int i;

    for (i=0; open_data_codes[i].s != NULL; i++) {
        if (!strcmp(s, open_data_codes[i].s)) {
            return open_data_codes[i].c;
        }
    }

    return 0;
}

void open_data (GtkAction *action)
{
    if (!dataset_locked()) {
        int code = open_data_code(gtk_action_get_name(action));

        file_selector(code, FSEL_DATA_NONE, NULL);
    }
}

void file_save (windata_t *vwin, int ci)
{
    switch (ci) {
    case SAVE_OUTPUT:
    case SAVE_CONSOLE:
    case SAVE_SCRIPT:
    case SAVE_GP_CMDS:
    case SAVE_R_CMDS:
    case SAVE_OX_CMDS:
    case SAVE_OCTAVE_CMDS:
    case SAVE_PYTHON_CMDS:
    case SAVE_STATA_CMDS:
    case SAVE_JULIA_CODE:
    case SAVE_DYNARE_CODE:
    case SAVE_LPSOLVE_CODE:
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case SAVE_SPEC_FILE:
    case SAVE_X13_SPC:
    case SAVE_HELP_TEXT:
        file_selector(ci, FSEL_DATA_VWIN, vwin);
        break;
    case EXPORT_CSV:
    case EXPORT:
        data_export_selection_wrapper(ci);
        break;
    case SAVE_TEX:
    case SAVE_TEXT:
        file_selector(ci, FSEL_DATA_MISC, vwin->data);
        break;
    default:
        dummy_call();
    }
}

static int fsave_code (const gchar *s)
{
    if (!strcmp(s, "SaveDataAs"))
        return SAVE_DATA_AS;
    if (!strcmp(s, "ExportData"))
        return EXPORT;

    return SAVE_DATA;
}

void fsave_callback (GtkAction *action, gpointer p)
{
    const gchar *s = gtk_action_get_name(action);
    int ci = fsave_code(s);

    file_save(p, ci);
}

static int model_action_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    int ci = gretl_command_number(s);

    if (ci == 0) {
        /* look up "GUI special" */
        if (!strcmp(s, "PANEL_WLS"))
            ci = PANEL_WLS;
        else if (!strcmp(s, "PANEL_B"))
            ci = PANEL_B;
        else if (!strcmp(s, "ALAGSEL"))
            ci = ALAGSEL;
        else if (!strcmp(s, "VLAGSEL"))
            ci = VLAGSEL;
        else if (!strcmp(s, "blogit"))
            ci = LOGIT;
        else if (!strcmp(s, "ologit"))
            ci = OLOGIT;
        else if (!strcmp(s, "mlogit"))
            ci = MLOGIT;
        else if (!strcmp(s, "bprobit"))
            ci = PROBIT;
        else if (!strcmp(s, "oprobit"))
            ci = OPROBIT;
        else if (!strcmp(s, "reprobit"))
            ci = REPROBIT;
        else if (!strcmp(s, "iv-liml"))
            ci = IV_LIML;
        else if (!strcmp(s, "iv-gmm"))
            ci = IV_GMM;
        else if (!strcmp(s, "countmod"))
            ci = COUNTMOD;
        else if (!strcmp(s, "regls"))
            ci = REGLS;
        else if (!strcmp(s, "FE_LOGISTIC"))
            ci = FE_LOGISTIC;
    }

    return ci;
}

void fit_resid_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) data;
    int code = M_UHAT;

    if (!strcmp(s, "yhat")) {
        code = M_YHAT;
    } else if (!strcmp(s, "uhat")) {
        code = M_UHAT;
    } else if (!strcmp(s, "uhat2")) {
        code = M_UHAT2;
    } else if (!strcmp(s, "h")) {
        code = M_H;
    } else if (!strcmp(s, "ahat")) {
        code = M_AHAT;
    }

    save_fit_resid(vwin, code);
}

/* callback for adding a scalar from a model, also
   pressed into service for saving a model as bundle
*/

void model_stat_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = vwin->data;
    int code = M_ESS;

    if (!strcmp(s, "ess")) {
        code = M_ESS;
    } else if (!strcmp(s, "se")) {
        code = M_SIGMA;
    } else if (!strcmp(s, "rsq")) {
        code = M_RSQ;
    } else if (!strcmp(s, "trsq")) {
        code = M_TRSQ;
    } else if (!strcmp(s, "df")) {
        code = M_DF;
    } else if (!strcmp(s, "lnL")) {
        code = M_LNL;
    } else if (!strcmp(s, "AIC")) {
        code = M_AIC;
    } else if (!strcmp(s, "BIC")) {
        code = M_BIC;
    } else if (!strcmp(s, "HQC")) {
        code = M_HQC;
    } else if (!strcmp(s, "bundle")) {
        code = B_MODEL;
    }

    add_model_stat(pmod, code, vwin);
}

static int have_midas_data (void)
{
    int i, m, got_midas = 0;

    for (i=1; i<dataset->v; i++) {
        m = series_is_midas_anchor(dataset, i);
        if (m > 0 && i + m <= dataset->v) {
            int is_midas = 1;
            int j, p, p0 = m;

            for (j=i+1; j<i+m; j++) {
                p = series_get_midas_period(dataset, j);
                if (p != p0 - 1) {
                    is_midas = 0;
                    break;
                } else {
                    p0 = p;
                }
            }
            if (is_midas) {
                got_midas = 1;
                break;
            }
        }
    }

    if (!got_midas) {
        const gchar *title = N_("gretl: warning");
        GtkWidget *dialog, *hbox, *help;
        const gchar *msg;

        msg = N_("No MIDAS data were found in the current dataset");
        dialog = gtk_message_dialog_new(GTK_WINDOW(mdata->main),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_MESSAGE_WARNING,
                                        GTK_BUTTONS_CLOSE,
                                        "%s", _(msg));
        g_signal_connect_swapped(dialog, "response",
                                 G_CALLBACK(gtk_widget_destroy),
                                 dialog);
        hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
        help = context_help_button(hbox, MIDAS_LIST);
        g_signal_connect_swapped(help, "clicked",
                                 G_CALLBACK(gtk_widget_destroy),
                                 dialog);
        gtk_widget_show(help);
        gtk_window_set_title(GTK_WINDOW(dialog), _(title));
        gtk_widget_show(dialog);
    }

    return got_midas;
}

void model_callback (GtkAction *action, gpointer data)
{
    int code = model_action_code(action);

    if (code == MIDASREG && !have_midas_data()) {
        return;
    } else {
        modelspec_dialog(code);
    }
}

void model_genr_callback (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    edit_dialog(MODEL_GENR, _("gretl: add var"),
                _("Enter formula for new variable:"),
                "", do_model_genr, vwin,
                VARCLICK_INSERT_NAME,
                vwin_toplevel(vwin));
}

static int selector_callback_code (const gchar *s)
{
    int ci = gretl_command_number(s);

    if (ci > 0) return ci;

    if (!strcmp(s, "TSPlot"))
        return GR_PLOT;
    if (!strcmp(s, "ScatterPlot"))
        return GR_XY;
    if (!strcmp(s, "ImpulsePlot"))
        return GR_IMP;
    if (!strcmp(s, "FactorPlot"))
        return GR_DUMMY;
    if (!strcmp(s, "FrischPlot"))
        return GR_XYZ;
    if (!strcmp(s, "ThreeDPlot"))
        return GR_3D;
    if (!strcmp(s, "MultiXY"))
        return SCATTERS;
    if (!strcmp(s, "MultiTS"))
        return TSPLOTS;
    if (!strcmp(s, "GR_BOX"))
        return GR_BOX;
    if (!strcmp(s, "GR_FBOX"))
        return GR_FBOX;
    if (!strcmp(s, "GR_QQ"))
        return QQPLOT;
    if (!strcmp(s, "ALAGSEL"))
        return ALAGSEL;
    if (!strcmp(s, "VLAGSEL"))
        return VLAGSEL;
    if (!strcmp(s, "ConfEllipse"))
        return ELLIPSE;
    if (!strcmp(s, "VarOmit"))
        return VAROMIT;
    if (!strcmp(s, "loess"))
        return LOESS;
    if (!strcmp(s, "nadarwat"))
        return NADARWAT;
    if (!strcmp(s, "Fsummary"))
        return FSUMMARY;

    return 0;
}

void selector_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) data;
    int ci;

    ci = selector_callback_code(s);

    if (ci == COINT || ci == COINT2) {
        selection_dialog(ci, _("gretl: cointegration test"), NULL, do_coint);
    } else if (ci == VAR || ci == VECM) {
        selection_dialog(ci, (ci == VAR)? _("gretl: VAR") : _("gretl: VECM"),
                         NULL, do_vector_model);
    } else if (ci == VLAGSEL) {
        selection_dialog(ci, _("gretl: VAR lag selection"), NULL, do_vector_model);
    } else if (ci == ALAGSEL) {
        selection_dialog(ci, _("gretl: ARIMA lag selection"), NULL, do_model);
    } else if (ci == GR_XY || ci == GR_IMP || ci == GR_DUMMY ||
               ci == SCATTERS || ci == GR_3D || ci == GR_XYZ ||
               ci == GR_FBOX || ci == FSUMMARY) {
        int (*selfunc)() = NULL;

        switch (ci) {
        case GR_XY:
        case GR_IMP:
            selfunc = do_graph_from_selector;
            break;
        case GR_3D:
            selfunc = do_splot_from_selector;
            break;
        case GR_DUMMY:
            selfunc = do_dummy_graph;
            break;
        case GR_XYZ:
            selfunc = do_xyz_graph;
            break;
        case SCATTERS:
            selfunc = do_multi_plots;
            break;
        case GR_FBOX:
        case FSUMMARY:
            selfunc = do_factorized_command;
            break;
        default:
            return;
        }
        if (ci == FSUMMARY) {
            selection_dialog(ci, _("gretl: factorized statistics"), NULL, selfunc);
        } else {
            selection_dialog(ci, _("gretl: define graph"), NULL, selfunc);
        }
    } else if (ci == ADD || ci == OMIT) {
        simple_selection_for_viewer(ci, _("gretl: model tests"),
                                    do_add_omit, vwin);
    } else if (ci == VAROMIT) {
        simple_selection_for_viewer(ci, _("gretl: model tests"),
                                    do_VAR_omit, vwin);
    } else if (ci == COEFFSUM) {
        simple_selection_for_viewer(ci, _("gretl: model tests"),
                                    do_coeff_sum, vwin);
    } else if (ci == ELLIPSE) {
        simple_selection_for_viewer(ci, _("gretl: model tests"),
                                    do_confidence_region, vwin);
    } else if (ci == GR_PLOT) {
        simple_selection(ci, _("gretl: define graph"), do_graph_from_selector, NULL);
    } else if (ci == TSPLOTS) {
        simple_selection(ci, _("gretl: define graph"), do_multi_plots, NULL);
    } else if (ci == GR_BOX) {
        simple_selection(ci, _("gretl: define graph"), do_regular_boxplot, NULL);
    } else if (ci == QQPLOT) {
        simple_selection(ci, _("gretl: define graph"), do_qq_from_selector, NULL);
    } else if (ci == LOESS || ci == NADARWAT) {
        const char *trstr = (ci == LOESS)? N_("Loess") : N_("Nadaraya-Watson");
        gchar *title;

        title = g_strdup_printf("gretl: %s", _(trstr));
        selection_dialog(ci, title, NULL, do_nonparam_model);
        g_free(title);
    } else {
        errbox("selector_callback: code was not recognized");
    }
}

static int gretl_callback_code (const gchar *s)
{
    if (!strcmp(s, "SampleRestrict"))
        return SMPLBOOL;
    if (!strcmp(s, "GENR"))
        return GENR;
    if (!strcmp(s, "VSETMISS"))
        return VSETMISS;
    if (!strcmp(s, "GSETMISS"))
        return GSETMISS;
    if (!strcmp(s, "gmm"))
        return GMM;
    if (!strcmp(s, "mle"))
        return MLE;
    if (!strcmp(s, "nls"))
        return NLS;
    if (!strcmp(s, "system"))
        return SYSTEM;
    if (!strcmp(s, "kalman"))
        return KALMAN;
    if (!strcmp(s, "restrict"))
        return RESTRICT;
    if (!strcmp(s, "MINIBUF"))
        return MINIBUF;
    return 0;
}

/* callback for menu items, where we want to prepare an "edit dialog"
   to handle the request */

void gretl_callback (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    GtkWidget *parent = NULL;
    const char *title = NULL;
    const char *query = NULL;
    const char *defstr = NULL;
    gchar *dynquery = NULL;
    void (*okfunc)() = NULL;
    Varclick click = VARCLICK_NONE;
    int ci;

    ci = gretl_callback_code(gtk_action_get_name(action));

    switch (ci) {
    case GENR:
        title = N_("gretl: add var");
        if (dataset->n > 5000) {
            query = N_("Enter formula for new variable");
        } else {
            query = N_("Enter formula for new variable\n"
                       "(or just a name, to enter data manually)");
        }
        okfunc = do_genr;
        click = VARCLICK_INSERT_NAME;
        break;
    case VSETMISS:
        title = N_("gretl: missing code");
        dynquery = g_strdup_printf(_("Enter value to be read as \"missing\"\n"
                                     "for the variable \"%s\""),
                                   dataset->varname[mdata->active_var]);
        okfunc = do_variable_setmiss;
        break;
    case GSETMISS:
        title = N_("gretl: missing code");
        query = N_("Enter value to be read as \"missing\":");
        okfunc = do_global_setmiss;
        break;
    case GMM:
        title = N_("gretl: GMM");
        query = N_("GMM: Specify function and orthogonality conditions:");
        okfunc = do_gmm_model;
        click = VARCLICK_INSERT_TEXT;
        data = NULL;
        break;
    case MLE:
        title = N_("gretl: maximum likelihood");
        query = N_("MLE: Specify function, and derivatives if possible:");
        okfunc = do_mle_model;
        click = VARCLICK_INSERT_TEXT;
        data = NULL;
        break;
    case NLS:
        title = N_("gretl: nonlinear least squares");
        query = N_("NLS: Specify function, and derivatives if possible:");
        okfunc = do_nls_model;
        click = VARCLICK_INSERT_TEXT;
        data = NULL;
        break;
    case SYSTEM:
        title = N_("gretl: simultaneous equations system");
        query = N_("Specify simultaneous equations:");
        data = NULL;
        okfunc = do_eqn_system;
        click = VARCLICK_INSERT_TEXT;
        break;
    case RESTRICT:
        title = N_("gretl: linear restrictions");
        query = N_("Specify restrictions:");
        okfunc = do_restrict;
        parent = vwin_toplevel(vwin);
        break;
    case MINIBUF:
        title = N_("gretl: command entry");
        query = N_("Type a command:");
        okfunc = do_minibuf;
        break;
    default:
        fprintf(stderr, "gretl_callback: unrecognized action '%s'\n",
               gtk_action_get_name(action));
        return;
    }

    if (dynquery != NULL) {
        edit_dialog(ci, _(title), dynquery, defstr, okfunc, data,
                    click, parent);
        g_free(dynquery);
    } else {
        edit_dialog(ci, _(title), _(query), defstr, okfunc, data,
                    click, parent);
    }
}

void menu_boxplot_callback (int varnum)
{
    const char *opts[] = {
        N_("Simple boxplot"),
        N_("With confidence interval for median"),
        N_("Factorized")
    };
    int ret;

    ret = radio_dialog(_("gretl: define graph"), NULL,
                       opts, 3, 0, 0, NULL);

    if (ret == 0) {
        do_boxplot_var(varnum, OPT_NONE);
    } else if (ret == 1) {
        do_boxplot_var(varnum, OPT_O);
    } else if (ret == 2) {
        selector_set_varnum(varnum);
        selection_dialog(GR_FBOX, _("gretl: define graph"),
                         NULL, do_factorized_command);
    }
}

void boxplot_callback (void)
{
    menu_boxplot_callback(mdata_active_var());
}

void revise_nl_model (MODEL *pmod, GtkWidget *parent)
{
    const char *title = NULL;
    const char *query = NULL;
    const char *defstr;
    void (*okfunc)() = NULL;

    switch (pmod->ci) {
    case GMM:
        title = N_("gretl: GMM");
        query = N_("GMM: Specify function and orthogonality conditions:");
        okfunc = do_gmm_model;
        break;
    case MLE:
        title = N_("gretl: maximum likelihood");
        query = N_("MLE: Specify function, and derivatives if possible:");
        okfunc = do_mle_model;
        break;
    case NLS:
        title = N_("gretl: nonlinear least squares");
        query = N_("NLS: Specify function, and derivatives if possible:");
        okfunc = do_nls_model;
        break;
    }

    defstr = gretl_model_get_data(pmod, "nlinfo");

    edit_dialog(pmod->ci, _(title), _(query), defstr, okfunc, pmod,
                VARCLICK_INSERT_TEXT, parent);
}

void revise_system_model (void *ptr, GtkWidget *parent)
{
    equation_system *sys = (equation_system *) ptr;

    edit_dialog(SYSTEM, _("gretl: simultaneous equations system"),
                _("Specify simultaneous equations:"), NULL,
                do_eqn_system, sys,
                VARCLICK_INSERT_TEXT,
                parent);
}

void genr_callback (void)
{
    const char *msg;

    if (dataset->n > 5000) {
        msg = N_("Enter formula for new variable");
    } else {
        msg = N_("Enter formula for new variable\n"
                 "(or just a name, to enter data manually)");
    }

    edit_dialog(GENR, _("gretl: add var"), _(msg),
                NULL, do_genr, NULL,
                VARCLICK_INSERT_NAME, NULL);
}

void minibuf_callback (void)
{
    edit_dialog(MINIBUF, _("gretl: command entry"),
                _("Type a command:"), NULL,
                do_minibuf, NULL,
                VARCLICK_NONE, NULL);
}

void newdata_callback (void)
{
    if (!dataset_locked()) {
        new_data_structure_dialog();
    }
}

void edit_gfn_callback (void)
{
    gchar *syspath = NULL;
    gchar *dotpath = NULL;
    gchar *startdir = NULL;
    int edit_sys_ok = 0;
    int edit_dot_ok = 0;
    int n_opts = 1;

    syspath = g_strdup_printf("%sfunctions", gretl_home());
    if (gretl_write_access(syspath) == 0) {
        edit_sys_ok = 1;
        n_opts++;
    }

#ifdef OS_OSX
    dotpath = g_strdup_printf("%sfunctions", gretl_app_support_dir());
    if (gretl_write_access(dotpath) == 0) {
        edit_dot_ok = 1;
        n_opts++;
    }
#else
    dotpath = g_strdup_printf("%sfunctions", gretl_dotdir());
    if (gretl_write_access(dotpath) == 0) {
        edit_dot_ok = 1;
        n_opts++;
    }
#endif

    if (n_opts > 1) {
        const char *opts[n_opts];
        int resp;

        if (edit_sys_ok && edit_dot_ok) {
            opts[0] = N_("system gfn directory");
            opts[1] = N_("personal gfn directory");
            opts[2] = N_("current working directory");
        } else if (edit_sys_ok) {
            opts[0] = N_("system gfn directory");
            opts[1] = N_("current working directory");
        } else if (edit_dot_ok) {
            opts[0] = N_("personal gfn directory");
            opts[1] = N_("current working directory");
        }

        resp = radio_dialog(_("Open .gfn file"), _("Start looking in:"),
                            opts, n_opts, 0, 0, NULL);
        if (resp < 0) {
            /* canceled */
            return;
        }

        if (edit_sys_ok && edit_dot_ok) {
            startdir = resp == 0 ? syspath : resp == 1 ? dotpath : NULL;
        } else if (edit_sys_ok) {
            startdir = resp == 0 ? syspath : NULL;
        } else if (edit_dot_ok) {
            startdir = resp == 0 ? dotpath : NULL;
        }
    }

    if (startdir != NULL) {
        file_selector_with_startdir(OPEN_GFN, startdir, NULL);
    } else {
        file_selector(OPEN_GFN, FSEL_DATA_NONE, NULL);
    }

    g_free(syspath);
    g_free(dotpath);
}

void install_pkg_callback (void)
{
    file_selector(INSTALL_PKG, FSEL_DATA_NONE, NULL);
}

void xcorrgm_callback (void)
{
    if (mdata_selection_count() == 2) {
        do_xcorrgm(NULL);
    } else {
        gchar *title = g_strdup_printf("gretl: %s", _("cross-correlogram"));

        simple_selection(XCORRGM, title, do_xcorrgm, NULL);
        g_free(title);
    }
}

void cond_number_callback (void)
{
    gretl_matrix *X = NULL;
    gretl_matrix *XX = NULL;
    gretl_matrix *chk = NULL;
    int *list = main_window_selection_as_list();
    int resp, err = 0;

    resp = yes_no_cancel_dialog(_("Collinearity check"),
                                _("Include a constant?"),
                                NULL);
    if (resp < 0) {
        /* canceled */
        return;
    }

    if (resp == GRETL_YES) {
        int *lplus = gretl_list_new(list[0] + 1);
        int i;

        lplus[1] = 0;
        for (i=2; i<=lplus[0]; i++) {
            lplus[i] = list[i-1];
        }
        free(list);
        list = lplus;
    }

    X = gretl_matrix_data_subset(list, dataset,
                                 dataset->t1,
                                 dataset->t2,
                                 M_MISSING_SKIP,
                                 &err);
    if (!err) {
        XX = gretl_matrix_XTX_new(X);
        if (XX == NULL) {
            err = E_ALLOC;
        } else {
            chk = gretl_matrix_copy(XX);
            err = gretl_invert_symmetric_matrix(XX);
            if (err == E_NOTPD && chk != NULL) {
                int r, sverr = 0;

                r = gretl_matrix_rank(chk, NADBL, &sverr);
                if (!sverr && r < XX->cols) {
                    infobox_printf(_("X'X is singular to machine precision:\n"
                                     "it is of dimension %d but rank %d."),
                                   XX->rows, r);
                    goto finish;
                }
            }
        }
    }

    if (err) {
        gui_errmsg(err);
    } else {
        gretl_matrix *(*bkwfunc) (const gretl_matrix *, gretl_array *,
                                  PRN *, int *);
        PRN *prn = NULL;

        bkwfunc = gui_get_plugin_function("bkw_matrix");
        if (bkwfunc != NULL && bufopen(&prn) == 0) {
            (*bkwfunc)(XX, NULL, prn, &err);
            if (err) {
                gui_errmsg(err);
            } else {
                view_buffer(prn, 78, 400, _("Collinearity check"), PRINT, NULL);
                /* record command? */
            }
        }
    }

 finish:

    gretl_matrix_free(X);
    gretl_matrix_free(XX);
    gretl_matrix_free(chk);
    free(list);
}

void map_save_callback (void)
{
    file_selector(WRITE_MAP, FSEL_DATA_NONE, NULL);
}

static int nist_verbosity (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "NistVVerbose"))
        return 2;
    else if (!strcmp(s, "NistVerbose"))
        return 1;
    else
        return 0;
}

void do_nistcheck (GtkAction *action)
{
    int (*run_nist_tests)(const char *, const char *, int);
    gchar *datadir = NULL;
    gchar *fname = NULL;

    run_nist_tests = gui_get_plugin_function("run_nist_tests");
    if (run_nist_tests == NULL) {
        return;
    }

    datadir = g_strdup_printf("%sdata%s", gretl_home(), SLASHSTR);
    fname = g_strdup_printf("%snist.out", gretl_dotdir());

    (*run_nist_tests)(datadir, fname, nist_verbosity(action));

    view_file(fname, 0, TMP_FILE, 78, 400, VIEW_CODEBOOK);

    g_free(datadir);
    g_free(fname);
}

static void mailer_help (void)
{
    show_gui_help(MAILHELP);
}

void send_attachment (const char *filename)
{
    int (*email_file) (const char *, GtkWindow *, void *);
    int err = 0;

    email_file = gui_get_plugin_function("email_file");

    if (email_file != NULL) {
        set_plugin_dialog_open(1);
        err = email_file(filename, GTK_WINDOW(mdata->main), mailer_help);
        set_plugin_dialog_open(0);
        if (err) {
            gui_errmsg(err);
        }
    }
}
