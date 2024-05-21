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

#include "gretl.h"
#include "var.h"
#include "johansen.h"
#include "varprint.h"
#include "forecast.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_func.h"
#include "system.h"
#include "matrix_extra.h"
#include "bootstrap.h"
#include "gretl_foreign.h"
#include "gretl_typemap.h"
#include "uservar.h"
#include "gretl_string_table.h"
#include "gretl_panel.h"
#include "csvdata.h"
#include "libset.h"

#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "model_table.h"
#include "series_view.h"
#include "session.h"
#include "textbuf.h"
#include "textutil.h"
#include "cmdstack.h"
#include "filelists.h"
#include "menustate.h"
#include "dlgutils.h"
#include "ssheet.h"
#include "datafiles.h"
#include "gpt_control.h"
#include "fileselect.h"
#include "toolbar.h"
#include "winstack.h"
#include "fnsave.h"
#include "datawiz.h"
#include "selector.h"
#include "guiprint.h"
#include "fncall.h"
#include "tabwin.h"
#include "join-gui.h"
#include "gretl_ipc.h"

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#endif

static void set_up_model_view_menu (windata_t *vwin);
static void add_system_menu_items (windata_t *vwin, int role);
static void add_x12_output_menu_item (windata_t *vwin);
static gint check_model_menu (GtkWidget *w, GdkEventButton *eb,
			      gpointer data);
static gint check_VAR_menu (GtkWidget *w, GdkEventButton *eb,
			    gpointer data);
static void model_copy_callback (GtkAction *action, gpointer p);
static int set_sample_from_model (void *ptr, int role);
static gboolean maybe_set_sample_from_model (windata_t *vwin);

static void close_model (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    gretl_viewer_destroy(vwin);
}

static int arma_by_x12a (const MODEL *pmod)
{
    int ret = 0;

    if (pmod->ci == ARMA) {
	int acode = gretl_model_get_int(pmod, "arma_flags");

	if (acode & ARMA_X12A) {
	    ret = 1;
	}
    }

    return ret;
}

static void model_output_save (GtkAction *action, gpointer p)
{
    copy_format_dialog((windata_t *) p, W_SAVE);
}

static gretlopt tex_eqn_opt;

static void set_tex_eqn_opt (GtkRadioAction *action)
{
    int v = gtk_radio_action_get_current_value(action);

    tex_eqn_opt = (v)? OPT_T : OPT_NONE;
}

gretlopt get_tex_eqn_opt (void)
{
    return tex_eqn_opt;
}

static int model_get_t1_t2 (void *ptr, int role, int *t1, int *t2)
{
    int err = 0;

    if (role == VIEW_MODEL) {
	MODEL *pmod = ptr;

	*t1 = pmod->smpl.t1;
	*t2 = pmod->smpl.t2;
    } else if (role == VAR || role == VECM) {
	GRETL_VAR *var = ptr;

	err = gretl_var_get_sample(var, t1, t2);
    } else if (role == SYSTEM) {
	equation_system *sys = ptr;

	err = gretl_system_get_sample(sys, t1, t2);
    }

    return err;
}

/* Called from menu in model window, but not necessarily
   a single-equation model */

static void model_revise_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    int ok = 1, t1 = 0, t2 = 0;
    int err = 0;

    if (vwin->role == VIEW_MODEL) {
	MODEL *pmod = vwin->data;

	err = model_sample_problem(pmod, dataset);
    }

    if (!err) {
	err = model_get_t1_t2(vwin->data, vwin->role, &t1, &t2);
    }

    if (!err && (t1 != dataset->t1 || t2 != dataset->t2)) {
	ok = maybe_set_sample_from_model(vwin);
    }

    if (ok) {
	selector_from_model(vwin);
    }
}

static void text_eqn_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    PRN *prn;
    int err;

    if (bufopen(&prn)) {
	return;
    }

    err = text_print_equation(pmod, dataset, OPT_NONE, prn);

    if (err) {
	gui_errmsg(err);
    } else {
	gchar *title = gretl_window_title(_("equation"));

	view_buffer_with_parent(vwin, prn, 78, 200, title, PRINT, NULL);
	g_free(title);
    }
}

static GtkActionEntry model_items[] = {
    { "File", NULL, N_("_File"), NULL, NULL, NULL },
    { "SaveAs", GTK_STOCK_SAVE_AS, N_("_Save as..."), NULL, NULL, G_CALLBACK(model_output_save) },
    { "SaveAsIcon", NULL, N_("Save to session as _icon"), NULL, NULL, G_CALLBACK(model_add_as_icon) },
    { "SaveAndClose", NULL, N_("Save as icon and cl_ose"), NULL, NULL, G_CALLBACK(model_add_as_icon) },
    { "Print", GTK_STOCK_PRINT, N_("_Print..."), NULL, NULL, G_CALLBACK(window_print) },
    { "TextEqn", NULL, N_("View as equation"), NULL, NULL, G_CALLBACK(text_eqn_callback) },
    { "Close", GTK_STOCK_CLOSE, N_("_Close"), NULL, NULL, G_CALLBACK(close_model) },
    { "Edit", NULL, N_("_Edit"), NULL, NULL, NULL },
    { "Copy", GTK_STOCK_COPY, N_("_Copy"), NULL, NULL, G_CALLBACK(model_copy_callback) },
    { "Revise", GTK_STOCK_EDIT, N_("_Modify model..."), NULL, NULL,
      G_CALLBACK(model_revise_callback) },
#if 0
    { "Restore", NULL, N_("_Restore model sample"), NULL, NULL, G_CALLBACK(model_sample_callback) },
#endif
    { "Tests", NULL, N_("_Tests"), NULL, NULL, NULL },
    { "Save", NULL, N_("_Save"), NULL, NULL, NULL },
    { "Graphs", NULL, N_("_Graphs"), NULL, NULL, NULL },
    { "ResidPlot", NULL, N_("_Residual plot"), NULL, NULL, NULL },
    { "FittedActualPlot", NULL, N_("_Fitted, actual plot"), NULL, NULL, NULL },
    { "Analysis", NULL, N_("_Analysis"), NULL, NULL, NULL },
    { "DisplayAFR", NULL, N_("_Display actual, fitted, residual"), NULL, NULL,
      G_CALLBACK(display_fit_resid) },
    { "Forecasts", NULL, N_("_Forecasts..."), NULL, NULL, G_CALLBACK(gui_do_forecast) },
    { "ConfIntervals", NULL, N_("_Confidence intervals for coefficients"), NULL, NULL,
      G_CALLBACK(do_coeff_intervals) },
    { "ConfEllipse", NULL, N_("Confidence _ellipse..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "Covariance", NULL, N_("Coefficient covariance _matrix"), NULL, NULL, G_CALLBACK(do_outcovmx) },
    { "Collinearity", NULL, N_("_Collinearity"), NULL, NULL, G_CALLBACK(do_collin) },
    { "Leverage", NULL, N_("_Influential observations"), NULL, NULL, G_CALLBACK(do_leverage) },
    { "ANOVA", NULL, N_("_ANOVA"), NULL, NULL, G_CALLBACK(do_anova) },
    { "Bootstrap", NULL, N_("_Bootstrap..."), NULL, NULL, G_CALLBACK(do_bootstrap) }
};

static GtkActionEntry model_test_items[] = {
    { "omit", NULL, N_("_Omit variables"), NULL, NULL, G_CALLBACK(selector_callback) },
    { "add", NULL, N_("_Add variables"), NULL, NULL, G_CALLBACK(selector_callback) },
    { "coeffsum", NULL, N_("_Sum of coefficients"), NULL, NULL, G_CALLBACK(selector_callback) },
    { "restrict", NULL, N_("_Linear restrictions"), NULL, NULL, G_CALLBACK(gretl_callback) },
    { "modtest:s", NULL, N_("Non-linearity (s_quares)"), NULL, NULL, G_CALLBACK(do_modtest) },
    { "modtest:l", NULL, N_("Non-linearity (_logs)"), NULL, NULL, G_CALLBACK(do_modtest) },
    { "reset", NULL, N_("_Ramsey's RESET"), NULL, NULL, G_CALLBACK(do_reset) },
    { "Hsk", NULL, N_("_Heteroskedasticity"), NULL, NULL, NULL },
    { "modtest:n", NULL, N_("_Normality of residual"), NULL, NULL, G_CALLBACK(do_resid_freq) },
    { "chow", NULL, N_("_Chow test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },
    { "modtest:a", NULL, N_("_Autocorrelation"), NULL, NULL, G_CALLBACK(do_autocorr) },
    { "dwpval", NULL, N_("_Durbin-Watson p-value"), NULL, NULL, G_CALLBACK(do_dwpval) },
    { "modtest:h", NULL, N_("A_RCH"), NULL, NULL, G_CALLBACK(do_arch) },
    { "bds", NULL, N_("Non-linearity (_BDS)"), NULL, NULL, G_CALLBACK(do_resid_freq) },
    { "qlrtest", NULL, N_("_QLR test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },
    { "cusum", NULL, N_("_CUSUM test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },
    { "cusum:r", NULL, N_("CUSUM_SQ test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },
    { "modtest:c", NULL, N_("_Common factor"), NULL, NULL, G_CALLBACK(do_modtest) },
    { "modtest:d", NULL, N_("_Cross-sectional dependence"), NULL, NULL, G_CALLBACK(do_modtest) },
    { "panspec", NULL, N_("_Panel specification"), NULL, NULL, G_CALLBACK(do_panel_tests) }
};

static GtkActionEntry base_hsk_items[] = {
    { "White", NULL, N_("White's test"), NULL, NULL, G_CALLBACK(do_modtest) },
    { "WhiteSquares", NULL, N_("White's test (squares only)"), NULL, NULL, G_CALLBACK(do_modtest) },
    { "BreuschPagan", NULL, "Breusch-Pagan", NULL, NULL, G_CALLBACK(do_modtest) },
    { "Koenker", NULL, "Koenker", NULL, NULL, G_CALLBACK(do_modtest) }
};

static GtkActionEntry panel_hsk_items[] = {
    { "White", NULL, N_("White's test"), NULL, NULL, G_CALLBACK(do_modtest) },
    { "Groupwise", NULL, N_("_groupwise"), NULL, NULL, G_CALLBACK(do_modtest) }
};

static GtkActionEntry ivreg_hsk_items[] = {
    { "White", NULL, N_("Pesaran-Taylor test"), NULL, NULL, G_CALLBACK(do_modtest) }
};

const gchar *model_tex_ui =
    "<ui>"
    "  <menubar>"
    "    <menu action='LaTeX'>"
    "      <menu action='TeXView'>"
    "        <menuitem action='TabView'/>"
    "        <menuitem action='EqnView'/>"
    "      </menu>"
    "      <menu action='TeXCopy'>"
    "        <menuitem action='TabCopy'/>"
    "        <menuitem action='EqnCopy'/>"
    "      </menu>"
    "      <menu action='TeXSave'>"
    "        <menuitem action='TabSave'/>"
    "        <menuitem action='EqnSave'/>"
    "      </menu>"
    "      <menu action='EqnOpts'>"
    "        <menuitem action='TeXstderrs'/>"
    "        <menuitem action='TeXtratios'/>"
    "      </menu>"
    "      <menuitem action='TabOpts'/>"
    "    </menu>"
    "  </menubar>"
    "</ui>";

const gchar *missing_tex_ui =
    "<ui>"
    "  <menubar>"
    "    <menu action='LaTeX'>"
    "      <menuitem action='notex'/>"
    "    </menu>"
    "  </menubar>"
    "</ui>";

static GtkActionEntry model_tex_items[] = {
    { "LaTeX",   NULL, N_("_LaTeX"), NULL, NULL, NULL },
    { "TeXView", NULL, N_("_View"), NULL, NULL, NULL },
    { "TabView", NULL, N_("_Tabular"), NULL, NULL, G_CALLBACK(model_tex_view) },
    { "EqnView", NULL, N_("_Equation"), NULL, NULL, G_CALLBACK(model_tex_view) },
    { "TeXCopy", NULL, N_("_Copy"), NULL, NULL, NULL },
    { "TabCopy", NULL, N_("_Tabular"), NULL, NULL, G_CALLBACK(model_tex_copy) },
    { "EqnCopy", NULL, N_("_Equation"), NULL, NULL, G_CALLBACK(model_tex_copy) },
    { "TeXSave", NULL, N_("_Save"), NULL, NULL, NULL },
    { "TabSave", NULL, N_("_Tabular"), NULL, NULL, G_CALLBACK(model_tex_save) },
    { "EqnSave", NULL, N_("_Equation"), NULL, NULL, G_CALLBACK(model_tex_save) },
    { "EqnOpts", NULL, N_("_Equation options"), NULL, NULL, NULL },
    { "TabOpts", NULL, N_("_Tabular options..."), NULL, NULL, G_CALLBACK(tex_format_dialog) }
};

static GtkRadioActionEntry tex_eqn_items[] = {
    { "TeXstderrs", NULL, N_("Show _standard errors"), NULL, NULL, 0 },
    { "TeXtratios", NULL, N_("Show _t-ratios"), NULL, NULL, 1 },
};

static GtkActionEntry missing_tex_items[] = {
    { "LaTeX", NULL, N_("_LaTeX"), NULL, NULL, NULL },
    { "notex", NULL, "No TeX", NULL, NULL, G_CALLBACK(dummy_call) }
};

static GtkActionEntry system_items[] = {
    { "File", NULL, N_("_File"), NULL, NULL, NULL },
    { "SaveAs", GTK_STOCK_SAVE_AS, N_("_Save as..."), NULL, NULL, G_CALLBACK(model_output_save) },
    { "SaveAsIcon", NULL, N_("Save to session as _icon"), NULL, NULL, G_CALLBACK(model_add_as_icon) },
    { "SaveAndClose", NULL, N_("Save as icon and cl_ose"), NULL, NULL, G_CALLBACK(model_add_as_icon) },
    { "Print", GTK_STOCK_PRINT, N_("_Print..."), NULL, NULL, G_CALLBACK(window_print) },
    { "Close", GTK_STOCK_CLOSE, N_("_Close"), NULL, NULL, G_CALLBACK(close_model) },
    { "Edit", NULL, N_("_Edit"), NULL, NULL, NULL },
    { "Copy", GTK_STOCK_COPY, N_("_Copy"), NULL, NULL, G_CALLBACK(model_copy_callback) },
    { "Revise", GTK_STOCK_EDIT, N_("_Revise specification..."), NULL, NULL,
      G_CALLBACK(model_revise_callback) },
    { "Save", NULL, N_("_Save"), NULL, NULL, NULL },
    { "Tests", NULL, N_("_Tests"), NULL, NULL, NULL },
    { "Graphs", NULL, N_("_Graphs"), NULL, NULL, NULL },
    { "Analysis", NULL, N_("_Analysis"), NULL, NULL, NULL },
    { "Forecasts", NULL, N_("_Forecasts"), NULL, NULL, NULL },
};

static gint n_system_items = G_N_ELEMENTS(system_items);

static GtkActionEntry sys_tex_items[] = {
    { "LaTeX",   NULL, N_("_LaTeX"), NULL, NULL, NULL },
    { "TeXView", NULL, N_("_View"),  NULL, NULL, G_CALLBACK(model_tex_view) },
    { "TeXCopy", NULL, N_("_Copy"),  NULL, NULL, G_CALLBACK(model_tex_copy) },
    { "TeXSave", NULL, N_("_Save"),  NULL, NULL, G_CALLBACK(model_tex_save) },
};

static const gchar *sys_ui =
    "<ui>"
    "  <menubar>"
    "    <menu action='File'>"
    "      <menuitem action='SaveAs'/>"
    "      <menuitem action='SaveAsIcon'/>"
    "      <menuitem action='SaveAndClose'/>"
    "      <menuitem action='Print'/>"
    "      <menuitem action='Close'/>"
    "    </menu>"
    "    <menu action='Edit'>"
    "      <menuitem action='Copy'/>"
    "      <menuitem action='Revise'/>"
    "    </menu>"
    "    <menu action='Tests'/>"
    "    <menu action='Save'/>"
    "    <menu action='Graphs'/>"
    "    <menu action='Analysis'>"
    "      <menu action='Forecasts'/>"
    "    </menu>"
    "  </menubar>"
    "</ui>";

static void model_copy_callback (GtkAction *action, gpointer p)
{
    copy_format_dialog((windata_t *) p, W_COPY);
}

void mark_dataset_as_modified (void)
{
    data_status |= MODIFIED_DATA;
    set_sample_label(dataset);

    if (session_file_is_open()) {
	mark_session_changed();
    }
}

static void gui_record_data_opening (const char *fname,
				     const int *list)
{
    const char *recname = (fname != NULL)? fname : datafile;

    if (data_status & BOOK_DATA) {
	/* don't print the (platform-dependent) full
	   path to a datafile packaged with gretl
	*/
	char basename[MAXSTR];

	gretl_basename(basename, recname, 0);
	lib_command_sprintf("open %s", basename);
    } else if (strchr(recname, ' ') != NULL) {
	lib_command_sprintf("open \"%s\"", recname);
    } else {
	lib_command_sprintf("open %s", recname);
    }

    if (list != NULL && list[0] == 3) {
	/* record spreadsheet parameters */
	char parm[32];

	if (list[1] > 1) {
	    sprintf(parm, " --sheet=%d", list[1]);
	    lib_command_strcat(parm);
	}
	if (list[2] > 0) {
	    sprintf(parm, " --coloffset=%d", list[2]);
	    lib_command_strcat(parm);
	}
	if (list[3] > 0) {
	    sprintf(parm, " --rowoffset=%d", list[3]);
	    lib_command_strcat(parm);
	}
    }

    record_command_verbatim();

    if (*datafile != '\0') {
	mkfilelist(FILE_LIST_DATA, datafile, 0);
    }
}

#define file_opened(f) (f & (DATAFILE_OPENED | OPENED_VIA_CLI | OPENED_VIA_SESSION))

static void real_register_data (int flags, const char *user_fname,
				const int *list)
{
    /* basic accounting */
    data_status |= HAVE_DATA;
    orig_vars = dataset->v;

    /* set appropriate data_status bits */
    if (file_opened(flags)) {
	if (!(data_status & IMPORT_DATA)) {
	    /* we opened a native data file */
	    if (has_system_prefix(datafile, DATA_SEARCH)) {
		data_status |= BOOK_DATA;
		data_status &= ~USER_DATA;
	    } else {
		data_status &= ~BOOK_DATA;
		data_status |= USER_DATA;
	    }
	    if (is_gzipped(datafile)) {
		data_status |= GZIPPED_DATA;
	    } else {
		data_status &= ~GZIPPED_DATA;
	    }
	    if (flags & OPENED_VIA_SESSION) {
		data_status |= SESSION_DATA;
	    } else {
		data_status &= ~SESSION_DATA;
	    }
	}
	if (flags & DATA_APPENDED) {
	    mark_dataset_as_modified();
	}
    } else {
	/* we modified the current dataset somehow */
	data_status |= GUI_DATA;
	mark_dataset_as_modified();
    }

    /* sync main window with datafile */
    if (mdata != NULL) {
	populate_varlist();
	set_sample_label(dataset);
	dataset_menubar_state(TRUE);
	session_menu_state(TRUE);
    }

    /* Record the opening of the data file in the GUI recent files
       list and command log; note that we don't do this if the file
       was opened via script or console, or if it was opened as a
       side effect of re-opening a saved session. And we can't do
       it if the data file was opened via the initial command line,
       and the gretl GUI is not yet built.
    */
    if (mdata != NULL && (flags & DATAFILE_OPENED)) {
	gui_record_data_opening(user_fname, list);
    }

    if (mdata != NULL) {
	if (!(flags & (OPENED_VIA_CLI | FOCUS_CONSOLE))) {
	    /* focus the data window */
	    gtk_widget_grab_focus(mdata->listbox);
	}
	/* invalidate "remove extra obs" menu item */
	drop_obs_state(FALSE);
    }
}

void register_data (int flags)
{
    real_register_data(flags, NULL, NULL);
}

void register_startup_data (const char *fname)
{
    real_register_data(DATAFILE_OPENED, fname, NULL);
}

static void maybe_offer_daily_options (GtkWidget *parent)
{
    gretlopt purge_opt = 0;
    PRN *prn = NULL;

    bufopen(&prn);

    if (prn != NULL) {
	int chk = analyse_daily_import(dataset, prn);
	const char *buf = NULL;
	int resp;

	if (chk > 0) {
	    buf = gretl_print_get_buffer(prn);
	}

	if (chk == 1) {
	    const char *opts[] = {
		N_("Leave the dataset as it is"),
		N_("Delete the weekend rows")
	    };

	    resp = radio_dialog("gretl", _(buf), opts, 2,
				1, DAILY_PURGE, parent);
	    if (resp == 1) {
		purge_opt = OPT_W;
	    }
	} else if (chk == 2) {
	    const char *opts[] = {
		N_("Leave the dataset as it is"),
		N_("Delete the weekend rows"),
		N_("Delete all blank rows")
	    };

	    resp = radio_dialog("gretl", _(buf), opts, 3,
				2, DAILY_PURGE, parent);
	    if (resp == 1) {
		purge_opt = OPT_W;
	    } else if (resp == 2) {
		purge_opt = OPT_A;
	    }
	} else if (chk == 3) {
	    const char *opts[] = {
		N_("Leave the dataset as it is"),
		N_("Delete all blank rows")
	    };

	    resp = radio_dialog("gretl", _(buf), opts, 2,
				1, DAILY_PURGE, parent);
	    if (resp == 1) {
		purge_opt = OPT_A;
	    }
	}

	gretl_print_destroy(prn);
    }

    if (purge_opt) {
	purge_opt |= OPT_T;
	bool_subsample(NULL, purge_opt, NULL);
    }
}

static void finalize_data_open (const char *fname, int ftype,
				int import, int append,
				const int *plist,
				GtkWidget *parent)
{
    if (import) {
	data_status |= IMPORT_DATA;
    }

    if (append) {
	register_data(DATA_APPENDED);
	return;
    }

    if (strstr(fname, CLIPTEMP_TXT) || strstr(fname, CLIPTEMP_GDT)) {
	real_register_data(DATA_PASTED, NULL, plist);
    } else {
	if (fname != datafile) {
	    strcpy(datafile, fname);
	}
	real_register_data(DATAFILE_OPENED, NULL, plist);
    }

    if (import && !dataset_is_time_series(dataset) &&
	!dataset_is_panel(dataset) && ftype != GRETL_MAP &&
	mdata != NULL) {
	int resp;

	resp = yes_no_dialog(_("gretl: open data"),
			     _("The imported data have been interpreted as undated\n"
			       "(cross-sectional).  Do you want to give the data a\n"
			       "time-series or panel interpretation?"),
			     parent);
	if (resp == GRETL_YES) {
	    data_structure_dialog();
	}
    } else if (import && dated_daily_data(dataset) &&
	       !dataset->markers) {
	maybe_offer_daily_options(parent);
    }
}

static int datafile_missing (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");
    int err = 0;

    if (fp == NULL) {
	delete_from_filelist(FILE_LIST_DATA, fname);
	file_read_errbox(fname);
	err = E_FOPEN;
    } else {
	fclose(fp);
    }

    return err;
}

/* below: get data of a sort that requires an import plugin */

int get_imported_data (char *fname, int ftype, int append,
		       gretlopt aopt)
{
    PRN *prn = NULL;
    int list[4] = {3, 0, 0, 0};
    int *plist = NULL;
    int (*ss_importer) (const char *, int *, char *, DATASET *,
			gretlopt, PRN *, GtkWidget *);
    int (*misc_importer) (const char *, DATASET *,
			  gretlopt, PRN *);
    gretlopt opt = OPT_G;
    int err = 0;

    if (datafile_missing(fname)) {
	return E_FOPEN;
    }

    import_na_init();
    ss_importer = NULL;
    misc_importer = NULL;

    if (ftype == GRETL_XLS) {
	ss_importer = gui_get_plugin_function("xls_get_data");
	plist = list;
    } else if (ftype == GRETL_XLSX) {
	ss_importer = gui_get_plugin_function("xlsx_get_data");
	plist = list;
    } else if (ftype == GRETL_GNUMERIC) {
	ss_importer = gui_get_plugin_function("gnumeric_get_data");
	plist = list;
    } else if (ftype == GRETL_ODS) {
	ss_importer = gui_get_plugin_function("ods_get_data");
	plist = list;
    } else if (ftype == GRETL_DTA) {
	misc_importer = gui_get_plugin_function("dta_get_data");
    } else if (ftype == GRETL_SAV) {
	misc_importer = gui_get_plugin_function("sav_get_data");
    } else if (ftype == GRETL_SAS) {
	misc_importer = gui_get_plugin_function("xport_get_data");
     } else if (ftype == GRETL_JMULTI) {
	misc_importer = gui_get_plugin_function("jmulti_get_data");
    } else if (ftype == GRETL_WF1) {
	misc_importer = gui_get_plugin_function("wf1_get_data");
    } else if (ftype == GRETL_MAP) {
	misc_importer = gui_get_plugin_function("map_get_data");
    } else {
	errbox(_("Unrecognized data type"));
	err = 1;
	goto bailout;
    }

    if (!err && ss_importer == NULL && misc_importer == NULL) {
	/* failed to open plugin */
        err = 1;
	goto bailout;
    }

    if (bufopen(&prn)) {
        err = 1;
	goto bailout;
    }

    if (append) {
	opt |= aopt;
    }

    /* call the actual importer function */
    if (SPREADSHEET_IMPORT(ftype)) {
	GtkWidget *parent = (mdata == NULL)? NULL : mdata->main;

	err = (*ss_importer)(fname, plist, NULL, dataset,
			     opt, prn, parent);
    } else {
	err = (*misc_importer)(fname, dataset, opt, prn);
    }

    if (err == -1) {
	fprintf(stderr, "data import canceled\n");
	err = E_CANCEL;
	goto bailout;
    } else {
	const char *buf = gretl_print_get_buffer(prn);

	if (err) {
	    if (buf != NULL && *buf != '\0') {
		errbox(buf);
	    } else {
		gui_errmsg(err);
	    }
	    if (err == E_FOPEN) {
		delete_from_filelist(FILE_LIST_DATA, fname);
	    }
	} else {
	    /* Do we want to be showing this? There's no error
	       but it might look sort of scary!
	    */
#if 0 /* masked out 2020-07-11 */
	    if (buf != NULL && *buf != '\0') {
		infobox(buf);
	    }
#endif
	    finalize_data_open(fname, ftype, 1, append, plist, NULL);
	}
    }

 bailout:

    gretl_print_destroy(prn);

    return err;
}

/* get "CSV" (or more generally, ASCII) data or GNU octave data:
   plugin is not required
*/

static int get_csv_data (char *fname, int ftype, int append,
			 gretlopt aopt)
{
    windata_t *vwin;
    PRN *prn;
    gchar *title = NULL;
    gretlopt opt = OPT_NONE;
    int err = 0;

    if (datafile_missing(fname)) {
	return E_FOPEN;
    }

    if (ftype == GRETL_CSV) {
        int resp = csv_open_dialog(fname);

	if (resp == GRETL_CANCEL) {
	    return 0;
	} else if (resp == 1) {
            /* all columns */
	    opt = OPT_A;
	}
    }

    if (bufopen(&prn)) {
	return 1;
    }

    if (append) {
	opt |= aopt;
    }

    if (ftype == GRETL_OCTAVE) {
	err = import_other(fname, ftype, dataset, opt, prn);
	title = g_strdup_printf(_("gretl: import %s data"), "Octave");
    } else {
	err = import_csv(fname, dataset, opt, prn);
	title = g_strdup_printf(_("gretl: import %s data"), "CSV");
    }

    /* show details regarding the import */
    vwin = view_buffer(prn, 78, 350, title, IMPORT, NULL);
    gtk_window_set_transient_for(GTK_WINDOW(vwin->main),
				 GTK_WINDOW(mdata->main));
    g_free(title);

    if (err) {
	delete_from_filelist(FILE_LIST_DATA, fname);
    } else {
	finalize_data_open(fname, ftype, 1, append, NULL,
			   vwin_toplevel(vwin));
    }

    return err;
}

static int get_native_data (char *fname, int ftype, int append,
			    gretlopt opt, windata_t *fwin)
{
    PRN *prn = NULL;
    char *buf = NULL;
    int err;

    if (bufopen(&prn)) {
	return 1;
    }

    if (ftype == GRETL_XML_DATA) {
	err = gretl_read_gdt(fname, dataset, opt | OPT_B, prn);
    } else {
	err = gretl_get_data(fname, dataset, opt, prn);
    }

    buf = gretl_print_steal_buffer(prn);
    gretl_print_destroy(prn);

    if (fwin != NULL) {
	/* close the files browser window that launched the query */
	gtk_widget_destroy(fwin->main);
    }

    if (err) {
	if (err == E_FOPEN) {
	    file_read_errbox(fname);
	} else if (buf != NULL && *buf != '\0') {
	    errbox(buf);
	} else {
	    gui_errmsg(err);
	}
	delete_from_filelist(FILE_LIST_DATA, fname);
    } else {
#if 0
	if (check_gretl_warning()) {
	    gui_warnmsg(0);
	}
#endif
	finalize_data_open(fname, ftype, 0, append, NULL, NULL);
	if (append) {
	    infobox(_("Data appended OK\n"));
	}
	fputs(buf, stderr);
    }

    free(buf);

    return err;
}

enum {
    SIMPLE_APPEND,
    FIXED_APPEND,
    USE_JOIN
};

static int select_append_type (void)
{
    const char *opts[] = {
	N_("simple append"),
	N_("append with fixed sample range"),
	N_("\"join\" (advanced)")
    };

    return radio_dialog(_("gretl: append data"), NULL,
		        opts, 3, 0, 0, NULL);
}

/* The gretl.c variable @tryfile will contain the name
   of the data file to be opened.
*/

gboolean do_open_data (windata_t *fwin, int code)
{
    int append = (code == APPEND_DATA);
    char *fname = get_tryfile();
    char tmp[MAXLEN];
    GretlFileType ftype;
    int append_type = SIMPLE_APPEND;
    gretlopt aopt = OPT_NONE;
    int err = 0;

    if (g_path_is_absolute(fname)) {
	strcpy(tmp, fname);
	ftype = data_file_type_from_name(tmp);
    } else {
	strcpy(tmp, fname);
	ftype = detect_filetype(tmp, OPT_P);
    }

    /* destroy the current data set, etc., unless we're
       explicitly appending */
    if (!append) {
	close_session(OPT_NONE); /* FIXME opt? */
    }

    if (append && (ftype == GRETL_CSV || ftype == GRETL_XML_DATA ||
		   ftype == GRETL_BINARY_DATA)) {
	append_type = select_append_type();
	if (append_type < 0) {
	    /* canceled */
	    return 0;
	} else if (append_type == FIXED_APPEND) {
	    aopt = OPT_X;
	}
    }

    if (append_type == USE_JOIN) {
	err = gui_join_data(tmp, ftype);
    } else if (ftype == GRETL_CSV || ftype == GRETL_OCTAVE) {
	err = get_csv_data(tmp, ftype, append, aopt);
    } else if (SPREADSHEET_IMPORT(ftype) || OTHER_IMPORT(ftype)) {
	err = get_imported_data(tmp, ftype, append, aopt);
    } else {
	err = get_native_data(tmp, ftype, append, aopt, fwin);
    }

    return !err;
}

static int give_dnd_options (GtkWidget *parent)
{
    const char *opts[] = {
	N_("Replace the current dataset (clears the current session)"),
	N_("Try appending to the current dataset"),
    };
    int resp;

    resp = radio_dialog(_("gretl: open data"), NULL, opts, 2,
			0, 0, parent);

    if (resp == GRETL_CANCEL) {
	return 0;
    } else {
	return resp == 0 ? OPEN_DATA : APPEND_DATA;
    }
}

static int give_cancel_option (GtkWidget *parent)
{
    const char *msg =
	N_("Opening a new data file will automatically\n"
	   "close the current one.  Any unsaved work\n"
	   "will be lost.  Proceed to open data file?");

    return yes_no_dialog(_("gretl: open data"), _(msg), parent);
}

/* give user choice of not opening selected datafile, if there's
   already a datafile open */

gboolean verify_open_data (windata_t *vwin, int action, gboolean dnd)
{
    if (dataset_locked()) {
	return FALSE;
    }

    if (data_status) {
	GtkWidget *parent = vwin_toplevel(vwin);

	if (dnd) {
	    action = give_dnd_options(parent);
	    if (action == 0) {
		return FALSE;
	    }
	} else {
	    int resp = give_cancel_option(parent);

	    if (resp != GRETL_YES) {
		return FALSE;
	    }
	}
    }

    return do_open_data(vwin, action);
}

/* give user choice of not opening session file, if there's already a
   datafile open */

gboolean verify_open_session (void)
{
    char *fname = get_tryfile();

    if (!gretl_is_pkzip_file(fname)) {
	/* not a zipped session file */
	return do_open_script(EDIT_HANSL);
    }

    if (data_status) {
	int resp =
	    yes_no_dialog(_("gretl: open session"),
			  _("Opening a new session file will automatically\n"
			    "close the current session.  Any unsaved work\n"
			    "will be lost.  Proceed to open session file?"),
			  NULL);

	if (resp != GRETL_YES) {
	    return FALSE;
	}
    }

    return do_open_session();
}

windata_t *console_window (int hsize, int vsize)
{
    windata_t *vwin;

    if (swallow) {
	vwin = gretl_viewer_new(CONSOLE, NULL, NULL);
    } else {
	gchar *title = NULL;
#ifdef GRETL_PID_FILE
	int seqno = gretl_sequence_number();

	if (seqno > 1) {
	    title = g_strdup_printf("%s (%d)", _("gretl console"), seqno);
	}
#endif
	if (title != NULL) {
	    vwin = gretl_viewer_new(CONSOLE, title, NULL);
	    g_free(title);
	} else {
	    vwin = gretl_viewer_new(CONSOLE, _("gretl console"), NULL);
	}
    }

    if (vwin == NULL) {
	return NULL;
    }

    vwin->flags |= VWIN_USE_FOOTER;
    vwin_add_viewbar(vwin, VIEWBAR_EDITABLE);
    create_console(vwin, hsize, vsize);
    text_table_setup(vwin->vbox, vwin->text);

    /* catch some special keystrokes */
    g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
		     G_CALLBACK(catch_viewer_key), vwin);

    gtk_widget_show(vwin->vbox);
    if (vwin->main != vwin->vbox) {
	gtk_widget_show(vwin->main);
    }

    g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
		     G_CALLBACK(text_popup_handler), vwin);

    return vwin;
}

windata_t *view_model (PRN *prn, MODEL *pmod, char *title)
{
    int hsize = MODEL_WIDTH;
    int vsize = MODEL_HEIGHT;
    windata_t *vwin;
    const char *buf;
    gchar *tmp;
    int width, nlines, tabbed;

#if 0
    fprintf(stderr, "view_model: pmod at %p, model %d\n",
	    (void *) pmod, pmod->ID);
#endif

    tabbed = use_tabbed_model_viewer();

    if (tabbed) {
	tmp = g_strdup_printf(_("model %d"), pmod->ID);
	vwin = viewer_tab_new(VIEW_MODEL, tmp, pmod);
	g_free(tmp);
    } else if (title != NULL) {
	vwin = gretl_viewer_new(VIEW_MODEL, title, pmod);
    } else {
	tmp = g_strdup_printf(_("gretl: model %d"), pmod->ID);
	vwin = gretl_viewer_new(VIEW_MODEL, tmp, pmod);
	g_free(tmp);
    }

    if (vwin == NULL) {
	return NULL;
    }

    /* Take responsibility for one reference to this model */
    gretl_object_ref(pmod, GRETL_OBJ_EQN);

    set_up_model_view_menu(vwin);

    gretl_print_get_size(prn, &width, &nlines);
    if (!tabbed && width > 0 && width + 2 < hsize) {
	hsize = width + 2;
    }

    create_text(vwin, hsize, vsize, nlines, FALSE);
    text_table_setup(vwin->vbox, vwin->text);

    /* insert and then free the model results buffer */
    buf = gretl_print_get_trimmed_buffer(prn);
    textview_set_text(vwin->text, buf);
    gretl_print_destroy(prn);

    /* record digits in force */
    widget_set_int(vwin->text, "digits", get_gretl_digits());

    /* sync number of model tests */
    vwin->n_model_tests = pmod->ntests;

    /* attach popup */
    g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
		     G_CALLBACK(text_popup_handler), vwin);

    if (window_is_tab(vwin)) {
	show_tabbed_viewer(vwin);
    } else {
	g_signal_connect(G_OBJECT(vwin->main), "key-press-event",
			 G_CALLBACK(catch_viewer_key), vwin);
	gtk_widget_show_all(vwin->main);
    }

    cursor_to_top(vwin);
    gtk_widget_grab_focus(vwin->text);

    return vwin;
}

static void mnl_probs_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    gretl_matrix *P = NULL;
    int err = 0;

    if (pmod == NULL) return;

    P = mn_logit_probabilities(pmod, pmod->t1, pmod->t2, dataset, &err);

    if (err) {
	gui_errmsg(err);
    } else {
	PRN *prn = NULL;

	if (bufopen(&prn)) {
	    gretl_matrix_free(P);
	} else {
	    gretl_matrix *yvals = gretl_model_get_data(pmod, "yvals");
	    int obslen = max_obs_marker_length(dataset);
	    int i, j, t = gretl_matrix_get_t1(P);
	    int n = gretl_vector_get_length(yvals);
	    double x;

	    pprintf(prn, _("\nEstimated outcome probabilities for %s\n\n"),
		    gretl_model_get_depvar_name(pmod, dataset));

	    /* case values */
	    bufspace(obslen, prn);
	    for (j=0; j<n; j++) {
		pprintf(prn, "%9g", yvals->val[j]);
	    }
	    pputc(prn, '\n');

	    /* format the matrix content nicely, prepending
	       observation strings */

	    for (i=0; i<P->rows; i++) {
		print_obs_marker(t, dataset, obslen, prn);
		for (j=0; j<P->cols; j++) {
		    x = gretl_matrix_get(P, i, j);
		    if (na(x)) {
			pprintf(prn, "%9s", " ");
		    } else {
			pprintf(prn, "%9.4f", x);
		    }
		}
		pputc(prn, '\n');
		t++;
	    }

	    view_buffer(prn, 78, 400, _("Outcome probabilities"), PRINT, NULL);
	    gretl_matrix_free(P);
	}
    }
}

static void add_multinomial_probs_item (windata_t *vwin)
{
    const gchar *mpath = "/menubar/Analysis";
    GtkActionEntry entry;

    action_entry_init(&entry);
    entry.name = "mnlprobs";
    entry.label = _("Outcome probabilities");
    entry.callback = G_CALLBACK(mnl_probs_callback);
    vwin_menu_add_item(vwin, mpath, &entry);
}

static void add_odds_ratios_item (windata_t *vwin)
{
    const gchar *mpath = "/menubar/Analysis";
    GtkActionEntry entry;

    action_entry_init(&entry);
    entry.name = "OddsRatios";
    entry.label = _("Odds ratios");
    entry.callback = G_CALLBACK(do_coeff_intervals);
    vwin_menu_add_item(vwin, mpath, &entry);
}

static int dw_pval_ok (const MODEL *pmod)
{
    if (na(pmod->dw)) {
	return 0;
    } else if (pmod->ci == OLS) {
	return 1;
    } else if (pmod->ci == PANEL) {
	return panel_DW_pval_ok(pmod);
    } else {
	return 0;
    }
}

static void get_ci_and_opt (const gchar *s, int *ci, gretlopt *opt)
{
    char c, word[9];

    sscanf(s, "%8[^:]:%c", word, &c);
    *ci = gretl_command_number(word);
    *opt = opt_from_flag((unsigned char) c);
}

static void set_tests_menu_state (GtkUIManager *ui, const MODEL *pmod)
{
    gretlopt opt;
    char path[128];
    const gchar *s;
    int i, n, ci;

    if (pmod->ci == MPOLS) {
	/* can we relax this? */
	flip(ui, "/menubar/Tests", FALSE);
	return;
    }

    n = G_N_ELEMENTS(model_test_items);

    for (i=0; i<n; i++) {
	int skip = 0;

	opt = OPT_NONE;
	s = model_test_items[i].name;
	if (strchr(s, ':')) {
	    get_ci_and_opt(s, &ci, &opt);
	} else if (!strcmp(s, "dwpval")) {
	    sprintf(path, "/menubar/Tests/%s", s);
	    flip(ui, path, dw_pval_ok(pmod));
	    continue;
	} else if (!strcmp(s, "Hsk")) {
	    ci = MODTEST;
	    if (pmod->ci == PANEL && (pmod->opt & OPT_U)) {
		/* random effects: not supported (FIXME?) */
		skip = 1;
	    } else {
		opt = (dataset_is_panel(dataset))? OPT_P : OPT_W;
	    }
	} else {
	    ci = gretl_command_number(s);
	}
	sprintf(path, "/menubar/Tests/%s", s);
	if (skip) {
	    flip(ui, path, FALSE);
	} else {
	    flip(ui, path, model_test_ok(ci, opt, pmod, dataset));
	}
    }

    if (pmod->ci == GARCH) {
	flip(ui, "/menubar/Tests/Hsk", FALSE);
    } else if (pmod->ci == PROBIT && (pmod->opt & OPT_E)) {
	/* random effects probit */
	flip(ui, "/menubar/Tests/modtest:n", FALSE);
    } else if (gretl_is_between_model(pmod)) {
	flip(ui, "/menubar/Tests/coeffsum", FALSE);
	flip(ui, "/menubar/Tests/Hsk", FALSE);
    }
}

static void set_analysis_menu_state (windata_t *vwin, const MODEL *pmod)
{
    GtkUIManager *ui = vwin->ui;
    gboolean s;

    if (pmod->ci == MLE || pmod->ci == GMM || pmod->ci == BIPROBIT) {
	/* can we relax some of this later? */
	flip(ui, "/menubar/Analysis/DisplayAFR", FALSE);
	flip(ui, "/menubar/Analysis/Forecasts", FALSE);
    } else if (pmod->ci == LOGIT) {
	if (gretl_model_get_int(pmod, "multinom")) {
	    /* relax this? */
	    flip(ui, "/menubar/Analysis/Forecasts", FALSE);
	    add_multinomial_probs_item(vwin);
	} else if (gretl_model_get_int(pmod, "binary")) {
	    add_odds_ratios_item(vwin);
	}
    } else if (pmod->ci == PROBIT && (pmod->opt & OPT_E)) {
	/* random effects probit */
	flip(ui, "/menubar/Analysis/Forecasts", FALSE);
    }

    if (pmod->ncoeff == 1) {
	flip(ui, "/menubar/Analysis/ConfEllipse", FALSE);
    }

    if (pmod->ci == DPANEL || (pmod->ci == PANEL && !(pmod->opt & OPT_P))) {
	flip(ui, "/menubar/Analysis/Forecasts", FALSE);
    }

    if (pmod->ci != OLS || !pmod->ifc || na(pmod->ess) || na(pmod->tss)) {
	flip(ui, "/menubar/Analysis/ANOVA", FALSE);
    }

    if (!bootstrap_ok(pmod->ci)) {
	flip(ui, "/menubar/Analysis/Bootstrap", FALSE);
    }

    if (gretl_model_get_int(pmod, "null-model")) {
	flip(ui, "/menubar/Analysis/ConfIntervals", FALSE);
	flip(ui, "/menubar/Analysis/Covariance", FALSE);
    }

    s = model_test_ok(VIF, OPT_NONE, pmod, dataset) ||
	pmod->vcv != NULL;
    flip(ui, "/menubar/Analysis/Collinearity", s);

    if (!model_test_ok(LEVERAGE, OPT_NONE, pmod, dataset)) {
	flip(ui, "/menubar/Analysis/Leverage", FALSE);
    }
}

static void arma_x12_menu_mod (windata_t *vwin)
{
    flip(vwin->ui, "/menubar/Analysis/Covariance", FALSE);
    add_x12_output_menu_item(vwin);
}

static void rq_coeff_intervals_mod (windata_t *vwin)
{
    flip(vwin->ui, "/menubar/Analysis/ConfIntervals", FALSE);
}

static void midas_plot_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    gretl_matrix *C = NULL;
    int err;

    if (pmod == NULL) return;

    C = gretl_model_get_data(pmod, "midas_coeffs");

    if (C != NULL) {
	const char *literal =
	    "{set title 'MIDAS coefficients';"
	    " set xlabel 'high-frequency lag';"
	    " set ylabel '';}";

	err = matrix_plot(C, NULL, literal,  OPT_P | OPT_S | OPT_G);
	gui_graph_handler(err);
    }
}

static void arma_spectrum_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int err;

    err = arma_spectrum_plot(pmod, dataset, OPT_NONE);

    gui_graph_handler(err);
}

#define intervals_model(m) (m->ci == LAD && \
			    gretl_model_get_data(m, "coeff_intervals"))

static void adjust_model_menu_state (windata_t *vwin, const MODEL *pmod)
{
    /* disallow saving an already-saved model */
    if (is_session_model((void *) pmod)) {
	set_model_save_state(vwin, FALSE);
    }

    if (RQ_SPECIAL_MODEL(pmod)) {
	/* can we relax this later? */
	flip(vwin->ui, "/menubar/Tests", FALSE);
	flip(vwin->ui, "/menubar/Save", FALSE);
	flip(vwin->ui, "/menubar/Analysis", FALSE);
	return;
    }

    set_tests_menu_state(vwin->ui, pmod);
    set_analysis_menu_state(vwin, pmod);

    if (intervals_model(pmod)) {
	rq_coeff_intervals_mod(vwin);
    }

    if (pmod->ci == ARMA && arma_by_x12a(pmod)) {
	arma_x12_menu_mod(vwin);
    }

    if (pmod->ci == MLE || pmod->ci == GMM || pmod->ci == BIPROBIT) {
	/* can we relax this later? */
	flip(vwin->ui, "/menubar/Graphs", FALSE);
    } else if (gretl_is_between_model(pmod) &&
	       pmod->dataset == NULL) {
	flip(vwin->ui, "/menubar/Graphs", FALSE);
    }

    if (pmod->ci == GMM || pmod->ci == BIPROBIT) {
	/* FIXME? */
	flip(vwin->ui, "/menubar/Save", FALSE);
    }
}

static GtkActionEntry model_data_base_items[] = {
    { "bundle", NULL, N_("_Model as bundle"), NULL, NULL,
      G_CALLBACK(model_stat_callback) },
    { "yhat", NULL, N_("_Fitted values"), NULL, NULL,
      G_CALLBACK(fit_resid_callback) },
    { "uhat", NULL, N_("_Residuals"), NULL, NULL,
      G_CALLBACK(fit_resid_callback) },
    { "uhat2", NULL, N_("_Squared residuals"), NULL, NULL,
      G_CALLBACK(fit_resid_callback) }
};

static GtkActionEntry ess_items[] = {
    { "ess", NULL, N_("_Error sum of squares"), NULL, NULL,
      G_CALLBACK(model_stat_callback) },
    { "se", NULL, N_("_Standard error of the regression"), NULL, NULL,
      G_CALLBACK(model_stat_callback) }
};

static GtkActionEntry r_squared_items[] = {
    { "rsq", NULL, N_("_R-squared"), NULL, NULL, G_CALLBACK(model_stat_callback) },
    { "trsq", NULL, N_("_T*R-squared"), NULL, NULL, G_CALLBACK(model_stat_callback) }
};

static GtkActionEntry lnl_data_items[] = {
    { "lnL", NULL, N_("_Log likelihood"), NULL, NULL,
      G_CALLBACK(model_stat_callback) }
};

static GtkActionEntry criteria_items[] = {
    { "AIC", NULL, N_("_Akaike Information Criterion"), NULL, NULL,
      G_CALLBACK(model_stat_callback) },
    { "BIC", NULL, N_("_Bayesian Information Criterion"), NULL, NULL,
      G_CALLBACK(model_stat_callback) },
    { "HQC", NULL, N_("_Hannan-Quinn Information Criterion"), NULL, NULL,
      G_CALLBACK(model_stat_callback) }
};

static GtkActionEntry garch_data_items[] = {
    { "h", NULL, N_("_Predicted error variance"), NULL, NULL,
      G_CALLBACK(fit_resid_callback)
    }
};

static GtkActionEntry fixed_effects_data_items[] = {
    { "ahat", NULL, N_("Per-unit _constants"), NULL, NULL,
      G_CALLBACK(fit_resid_callback)
    }
};

static GtkActionEntry random_effects_data_items[] = {
    { "ahat", NULL, N_("Individual effects"), NULL, NULL,
      G_CALLBACK(fit_resid_callback)
    }
};

static GtkActionEntry define_var_items[] = {
    /* Under Save; Sep wanted */
    { "NewVar", NULL, N_("Define _new variable..."), NULL, NULL,
      G_CALLBACK(model_genr_callback) }
};

static int criteria_available (const MODEL *pmod)
{
    int i;

    for (i=0; i<C_MAX; i++) {
	if (na(pmod->criterion[i])) {
	    return 0;
	}
    }

    return 1;
}

static void add_model_dataset_items (windata_t *vwin)
{
    const gchar *path = "/menubar/Save";
    MODEL *pmod = vwin->data;

    vwin_menu_add_items(vwin, path, model_data_base_items,
			G_N_ELEMENTS(model_data_base_items));

    if (gretl_model_get_data(pmod, "ahat") != NULL) {
	if (pmod->opt & OPT_U) {
	    vwin_menu_add_items(vwin, path, random_effects_data_items,
				G_N_ELEMENTS(random_effects_data_items));
	} else {
	    vwin_menu_add_items(vwin, path, fixed_effects_data_items,
				G_N_ELEMENTS(fixed_effects_data_items));
	}
    }

    if (pmod->ci != GARCH && !(pmod->ci == LOGIT && (pmod->opt & OPT_M))) {
	vwin_menu_add_items(vwin, path, ess_items,
			    G_N_ELEMENTS(ess_items));
    }

    if (!ML_ESTIMATOR(pmod->ci) && pmod->ci != LAD && !na(pmod->rsq)) {
	vwin_menu_add_items(vwin, path, r_squared_items,
			    G_N_ELEMENTS(r_squared_items));
    }

    if (!na(pmod->lnL)) {
	vwin_menu_add_items(vwin, path, lnl_data_items,
			    G_N_ELEMENTS(lnl_data_items));
    }

    if (criteria_available(pmod)) {
	vwin_menu_add_items(vwin, path, criteria_items,
			    G_N_ELEMENTS(criteria_items));
    }

    if (pmod->ci == GARCH) {
	vwin_menu_add_items(vwin, path, garch_data_items,
			    G_N_ELEMENTS(garch_data_items));
    }

    vwin_menu_add_separator(vwin, path);

    vwin_menu_add_items(vwin, path, define_var_items,
			G_N_ELEMENTS(define_var_items));
}

static void add_model_tex_items (windata_t *vwin)
{
    MODEL *pmod = (MODEL *) vwin->data;
    int eqn_ok = command_ok_for_model(EQNPRINT, 0, pmod);
    GtkActionGroup *actions;
    GError *err = NULL;
    int imod = 0;

    gtk_ui_manager_add_ui_from_string(vwin->ui, model_tex_ui, -1, &err);

    if (err != NULL) {
	g_message("building LaTeX menu failed: %s", err->message);
	g_error_free(err);
	return;
    }

    actions = gtk_action_group_new("ModelTeX");
    gtk_action_group_set_translation_domain(actions, "gretl");
    gtk_action_group_add_actions(actions, model_tex_items,
				 G_N_ELEMENTS(model_tex_items),
				 vwin);
    gtk_action_group_add_radio_actions(actions, tex_eqn_items,
				       G_N_ELEMENTS(tex_eqn_items),
				       (get_tex_eqn_opt() == OPT_T),
				       G_CALLBACK(set_tex_eqn_opt),
				       vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    if (intervals_model(pmod)) {
	eqn_ok = 0;
	imod = 1;
    }

    if (!eqn_ok || pmod->errcode) {
	flip(vwin->ui, "/menubar/LaTeX/TeXView/EqnView", FALSE);
	flip(vwin->ui, "/menubar/LaTeX/TeXCopy/EqnCopy", FALSE);
	flip(vwin->ui, "/menubar/LaTeX/TeXSave/EqnSave", FALSE);
	flip(vwin->ui, "/menubar/LaTeX/EqnOpts", FALSE);
    }

    if (imod) {
	flip(vwin->ui, "/menubar/LaTeX/TabOpts", FALSE);
    }
}

/* dummy placeholder, for when TeX is not supported */

static void add_missing_tex_items (windata_t *vwin)
{
    GtkActionGroup *actions;
    GError *err = NULL;

    gtk_ui_manager_add_ui_from_string(vwin->ui, missing_tex_ui, -1, &err);
    if (err != NULL) {
	g_message("building menus failed: %s", err->message);
	g_error_free(err);
	return;
    }

    actions = gtk_action_group_new("MissingTeX");
    gtk_action_group_add_actions(actions, missing_tex_items,
				 G_N_ELEMENTS(missing_tex_items),
				 vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    flip(vwin->ui, "/menubar/LaTeX", FALSE);
}

#define VNAMELEN2 (2*VNAMELEN)

static void add_vars_to_plot_menu (windata_t *vwin)
{
    GtkActionEntry entry;
    const gchar *mpath[] = {
	"/menubar/Graphs/ResidPlot",
	"/menubar/Graphs/FittedActualPlot"
    };
    MODEL *pmod = vwin->data;
    char tmp[VNAMELEN2], aname[VNAMELEN];
    gchar *alabel;
    int *xlist;
    int v1, v2;
    int i, j;

    action_entry_init(&entry);

    /* non-time series: residual boxplot */
    if (!dataset_is_time_series(dataset)) {
	entry.label = _("_Boxplot");
	entry.name = "r:box";
	entry.callback = G_CALLBACK(resid_plot);
	vwin_menu_add_item(vwin, mpath[0], &entry);
    }

    xlist = gretl_model_get_x_list(pmod);

    for (i=0; i<2; i++) {
	/* plot against time/obs number */
	if (dataset_is_time_series(dataset)) {
	    entry.label = _("_Against time");
	} else {
	    entry.label = _("By _observation number");
	}
	entry.name = (i == 0)? "r:byobs" : "f:byobs";
	entry.callback = (i == 0)? G_CALLBACK(resid_plot) :
	    G_CALLBACK(fit_actual_plot);
	vwin_menu_add_item(vwin, mpath[i], &entry);

	if (pmod->ci == NLS ||
	    pmod->ci == MLE ||
	    pmod->ci == GMM ||
	    pmod->ci == PANEL) {
	    continue;
	}

	/* if doing resid plot, put dependent var in menu */
	if (i == 0) {
	    v1 = gretl_model_get_depvar(pmod);
	    if (v1 > 0) {
		sprintf(aname, "r:xvar %d", v1); /* FIXME */
		double_underscores(tmp, dataset->varname[v1]);
		alabel = g_strdup_printf(_("_Against %s"), tmp);
		entry.name = aname;
		entry.label = alabel;
		entry.callback = G_CALLBACK(resid_plot);
		vwin_menu_add_item(vwin, mpath[0], &entry);
		g_free(alabel);
	    }
	}

	if (xlist != NULL) {
	    /* put the independent vars on the menu list */
	    for (j=1; j<=xlist[0]; j++) {
		v1 = xlist[j];
		if (v1 == 0) {
		    continue;
		}
		if (!strcmp(dataset->varname[v1], "time")) {
		    continue;
		}
		sprintf(aname, "%c:xvar %d", (i == 0)? 'r' : 'f', v1);
		double_underscores(tmp, dataset->varname[v1]);
		alabel = g_strdup_printf(_("_Against %s"), tmp);
		entry.name = aname;
		entry.label = alabel;
		entry.callback = (i == 0)? G_CALLBACK(resid_plot) :
		    G_CALLBACK(fit_actual_plot);
		vwin_menu_add_item(vwin, mpath[i], &entry);
		g_free(alabel);
	    }
	}

	if (i == 1) {
	    /* fitted values: offer Theil-type scatterplot */
	    entry.name = "f:theil";
	    entry.label = _("Actual vs. Fitted");
	    entry.callback = G_CALLBACK(fit_actual_plot);
	    vwin_menu_add_item(vwin, mpath[i], &entry);
	}
    }

    /* time series models: residual correlogram, spectrum */
    if (dataset_is_time_series(dataset)) {
	vwin_menu_add_separator(vwin, "/menubar/Graphs");
	entry.name = "Correlogram";
	entry.label = _("Residual _correlogram");
	entry.callback = G_CALLBACK(residual_correlogram_callback);
	vwin_menu_add_item(vwin, "/menubar/Graphs", &entry);
	entry.name = "Spectrum";
	entry.label = _("Residual _periodogram");
	entry.callback = G_CALLBACK(residual_periodogram_callback);
	vwin_menu_add_item(vwin, "/menubar/Graphs", &entry);
	if (pmod->ci == ARMA && gretl_model_get_data(pmod, "ainfo") != NULL) {
	    entry.name = "ARMAspectrum";
	    entry.label = _("_Spectrum vs sample periodogram");
	    entry.callback = G_CALLBACK(arma_spectrum_callback);
	    vwin_menu_add_item(vwin, "/menubar/Graphs", &entry);
	} else if (gretl_model_get_data(pmod, "midas_coeffs") != NULL) {
	    entry.name = "MIDAScoeffs";
	    entry.label = _("_MIDAS coefficients");
	    entry.callback = G_CALLBACK(midas_plot_callback);
	    vwin_menu_add_item(vwin, "/menubar/Graphs", &entry);
	}
    } else {
	vwin_menu_add_separator(vwin, "/menubar/Graphs");
    }

    /* residual Q-Q plot */
    entry.name = "QQPlot";
    entry.label = _("Residual _Q-Q plot");
    entry.callback = G_CALLBACK(residual_qq_plot);
    vwin_menu_add_item(vwin, "/menubar/Graphs", &entry);

    /* 3-D fitted versus actual plot? */
    if (xlist != NULL) {
	int parnames = pmod->dataset != NULL;

	v1 = v2 = -1;
	if (pmod->ifc && xlist[0] == 3) {
	    v1 = parnames ? 1 : xlist[2];
	    v2 = parnames ? 2 : xlist[3];
	} else if (!pmod->ifc && xlist[0] == 2) {
	    v1 = parnames ? 0 : xlist[1];
	    v2 = parnames ? 1 : xlist[2];
	}
	if (v1 >= 0 && v2 >= 0) {
	    char tmp2[VNAMELEN2];

	    if (parnames) {
		char pname[VNAMELEN];

		gretl_model_get_param_name(pmod, pmod->dataset,
					   v1, pname);
		double_underscores(tmp, pname);
		gretl_model_get_param_name(pmod, pmod->dataset,
					   v2, pname);
		double_underscores(tmp2, pname);
	    } else {
		double_underscores(tmp, dataset->varname[v1]);
		double_underscores(tmp2, dataset->varname[v2]);
	    }
	    vwin_menu_add_separator(vwin, mpath[1]);
	    alabel = g_strdup_printf(_("_Against %s and %s"), tmp, tmp2);
	    entry.name = "splot";
	    entry.label = alabel;
	    entry.callback = G_CALLBACK(fit_actual_splot);
	    vwin_menu_add_item(vwin, mpath[1], &entry);
	    g_free(alabel);
	}
    }

    free(xlist);
}

static void plot_dummy_call (GtkRadioAction *action,
			     GtkRadioAction *current,
			     windata_t *vwin)
{
    vwin->active_var = gtk_radio_action_get_current_value(action);
}

static void radio_action_init (GtkRadioActionEntry *a)
{
    a->stock_id = NULL;
    a->accelerator = NULL;
    a->tooltip = NULL;
}

static void add_dummies_to_plot_menu (windata_t *vwin)
{
    GtkActionEntry item;
    GtkRadioActionEntry *items;
    MODEL *pmod = vwin->data;
    const gchar *gpath = "/menubar/Graphs/ResidPlot";
    const gchar *spath = "/menubar/Graphs/ResidPlot/Separation";
    char tmp[VNAMELEN2];
    int *dlist = NULL;
    int i, vi, ndums;

    /* make a list of dummy independent variables */
    for (i=2; i<=pmod->list[0]; i++) {
	vi = pmod->list[i];
	if (vi == LISTSEP) {
	    break;
	} else if (vi > 0 && vi < dataset->v &&
	    gretl_isdummy(dataset->t1, dataset->t2, dataset->Z[vi])) {
	    gretl_list_append_term(&dlist, vi);
	}
    }

    if (dlist == NULL) {
	return;
    }

    ndums = dlist[0];
    items = malloc((ndums + 1) * sizeof *items);
    if (items == NULL) {
	free(dlist);
	return;
    }

    /* add separator */
    vwin_menu_add_separator(vwin, gpath);

    /* add menu branch */
    action_entry_init(&item);
    item.name = "Separation";
    item.label = _("Separation");
    vwin_menu_add_menu(vwin, gpath, &item);

    /* configure "none" radio option */
    radio_action_init(&items[0]);
    items[0].name = "none";
    items[0].label = _("none");
    items[0].value = 0;

    /* put the dummy independent vars on the menu list */
    for (i=1; i<=dlist[0]; i++) {
	vi = dlist[i];
	radio_action_init(&items[i]);
	double_underscores(tmp, dataset->varname[vi]);
	items[i].name = g_strdup_printf("dum %d", vi);
	items[i].label = g_strdup_printf(_("By %s"), tmp);
	items[i].value = vi;
    }

    vwin_menu_add_radios(vwin, spath, items, ndums + 1, 0,
			 G_CALLBACK(plot_dummy_call));

    for (i=1; i<=dlist[0]; i++) {
	g_free((gchar *) items[i].name);
	g_free((gchar *) items[i].label);
    }

    free(items);
    free(dlist);
}

static void varnum_from_action (GtkAction *action, int *i)
{
    const gchar *s = gtk_action_get_name(action);

    sscanf(s, "%*s %d", i);
}

static void tau_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int v, err;

    varnum_from_action(action, &v);
    err = plot_tau_sequence(pmod, dataset, v);
    gui_graph_handler(err);
}

static void add_tau_plot_menu (windata_t *vwin)
{
    GtkActionEntry item;
    MODEL *pmod = vwin->data;
    char tmp[VNAMELEN2], aname[VNAMELEN];
    int i;

    action_entry_init(&item);
    item.name = "TauMenu";
    item.label = _("tau sequence");
    vwin_menu_add_menu(vwin, "/menubar/Graphs", &item);

    item.callback = G_CALLBACK(tau_plot_call);

    /* put the independent vars on the menu list */
    for (i=2; i<=pmod->list[0]; i++) {
	sprintf(aname, "tauseq %d", i - 2);
	double_underscores(tmp, dataset->varname[pmod->list[i]]);
	item.name = aname;
	item.label = tmp;
	vwin_menu_add_item(vwin, "/menubar/Graphs/TauMenu", &item);
    }

    /* and disable what we can't (yet) show */
    flip(vwin->ui, "/menubar/Graphs/ResidPlot", FALSE);
    flip(vwin->ui, "/menubar/Graphs/FittedActualPlot", FALSE);
}

static void x12_output_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    char *fname;

    if (pmod == NULL) return;

    fname = gretl_model_get_data(pmod, "x12a_output");

    if (fname != NULL) {
	char *p = strrchr(fname, '.');

	if (p != NULL && strlen(p) == 7) {
	    gchar *tmp = g_strdup(fname);

	    sprintf(p, ".%d", pmod->ID);
	    gretl_rename(tmp, fname);
	    g_free(tmp);
	}
	view_file(fname, 0, 0, 78, 350, VIEW_FILE);
    }
}

static gchar *get_model_ui (void)
{
    gchar *ui = NULL;
    gchar *fname;
    int err;

    fname = g_strdup_printf("%sui%cgretlmodel.xml", gretl_home(), SLASH);
    err = gretl_file_get_contents(fname, &ui, NULL);
    g_free(fname);

    return err ? NULL : ui;
}

int latex_is_ok (void)
{
    static int latex_ok = -1;

    if (latex_ok == -1) {
	latex_ok = check_for_program(latex);
    }

    return latex_ok;
}

static void set_up_model_view_menu (windata_t *vwin)
{
    static gchar *model_ui;
    MODEL *pmod = (MODEL *) vwin->data;
    GtkActionGroup *actions;
    GtkWidget *toplevel;
    GError *err = NULL;

    if (model_ui == NULL) {
	model_ui = get_model_ui();
	if (model_ui == NULL) {
	    errbox("building menus failed");
	    return;
	}
    }

    actions = gtk_action_group_new("ModelActions");
    gtk_action_group_set_translation_domain(actions, "gretl");

    gtk_action_group_add_actions(actions, model_items,
				 G_N_ELEMENTS(model_items),
				 vwin);
    gtk_action_group_add_actions(actions, model_test_items,
				 G_N_ELEMENTS(model_test_items),
				 vwin);

    vwin->ui = gtk_ui_manager_new();
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    gtk_ui_manager_add_ui_from_string(vwin->ui, model_ui, -1, &err);
    if (err != NULL) {
	g_message("building menus failed: %s", err->message);
	g_error_free(err);
    }

    if (pmod->ci != MLE && pmod->ci != GMM && pmod->ci != BIPROBIT) {
	if (RQ_SPECIAL_MODEL(pmod)) {
	    add_tau_plot_menu(vwin);
	} else {
	    add_vars_to_plot_menu(vwin);
	}
	add_model_dataset_items(vwin);
    }

    /* heteroskedasticity tests: the permissible options vary
       depending on the nature of the model
    */

    if (dataset_is_panel(dataset)) {
	if (pmod->ci == OLS) {
	    vwin_menu_add_items(vwin, "/menubar/Tests/Hsk",
				panel_hsk_items,
				G_N_ELEMENTS(panel_hsk_items));
	} else if (pmod->ci == PANEL && (pmod->opt & OPT_F)) {
	    /* fixed effects */
	    vwin_menu_add_items(vwin, "/menubar/Tests/Hsk",
				panel_hsk_items + 1, 1);
	} else if (0 && pmod->ci == PANEL && (pmod->opt & OPT_U)) {
	    /* random effects: is this OK?? */
	    vwin_menu_add_items(vwin, "/menubar/Tests/Hsk",
				panel_hsk_items + 1, 1);
	}
    } else if (pmod->ci == IVREG) {
	vwin_menu_add_items(vwin, "/menubar/Tests/Hsk",
			    ivreg_hsk_items,
			    G_N_ELEMENTS(ivreg_hsk_items));
    } else if (model_test_ok(MODTEST, OPT_W, pmod, dataset)) {
	vwin_menu_add_items(vwin, "/menubar/Tests/Hsk",
			    base_hsk_items,
			    G_N_ELEMENTS(base_hsk_items));
	if (pmod->ncoeff == 1 || (pmod->ifc && pmod->ncoeff == 2)) {
	    flip(vwin->ui, "/menubar/Tests/Hsk/WhiteSquares", FALSE);
	}
    }

    maybe_add_packages_to_model_menus(vwin);

    if (latex_is_ok() && !pmod->errcode && !RQ_SPECIAL_MODEL(pmod)) {
	add_model_tex_items(vwin);
    } else {
	add_missing_tex_items(vwin);
    }

    if (!text_equation_ok(pmod)) {
	flip(vwin->ui, "/menubar/File/TextEqn", FALSE);
    }

    if (pmod->ci != ARMA && pmod->ci != GARCH &&
	pmod->ci != NLS && pmod->ci != MLE && pmod->ci != GMM &&
	pmod->ci != PANEL && pmod->ci != DPANEL &&
	pmod->ci != BIPROBIT) {
	add_dummies_to_plot_menu(vwin);
    }

    toplevel = vwin_toplevel(vwin);
    if (toplevel != NULL) {
	/* FIXME tabbed case? */
	gtk_window_add_accel_group(GTK_WINDOW(toplevel),
				   gtk_ui_manager_get_accel_group(vwin->ui));
    }

    vwin->mbar = gtk_ui_manager_get_widget(vwin->ui, "/menubar");

    /* disable some menu items if need be */
    adjust_model_menu_state(vwin, pmod);

    g_signal_connect(G_OBJECT(vwin->mbar), "button-press-event",
		     G_CALLBACK(check_model_menu), vwin);

    vwin_pack_toolbar(vwin);
}

enum {
    SYS_DATA_RESIDS,
    SYS_DATA_FITTED,
    SYS_DATA_SIGMA
};

static int sys_data_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "uhat")) {
	return SYS_DATA_RESIDS;
    } else if (!strcmp(s, "yhat")) {
	return SYS_DATA_FITTED;
    } else if (!strcmp(s, "sigma")) {
	return SYS_DATA_SIGMA;
    } else {
	return SYS_DATA_RESIDS;
    }
}

static void system_data_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const gretl_matrix *M = NULL;
    gchar *wtitle = NULL;
    PRN *prn;
    int code, k = 0;
    int err = 0;

    if (vwin->role == SYSTEM) {
	sys = (equation_system *) vwin->data;
    } else {
	var = (GRETL_VAR *) vwin->data;
    }

    if ((var == NULL && sys == NULL) || bufopen(&prn)) {
	return;
    }

    code = sys_data_code(action);

    if (code == SYS_DATA_SIGMA) {
	if (var != NULL) {
	    wtitle = g_strdup(_("gretl: VAR covariance matrix"));
	    err = gretl_VAR_print_sigma(var, prn);
	} else {
	    wtitle = g_strdup(_("gretl: system covariance matrix"));
	    err = system_print_sigma(sys, prn);
	}
    } else if (code == SYS_DATA_RESIDS || code == SYS_DATA_FITTED) {
	const char *titles[] = {
	    N_("System residuals"),
	    N_("System fitted values")
	};
	const char *title;
	const char **heads = NULL;

	if (var != NULL) {
	    /* fitted values matrix not currently available */
	    M = (code == SYS_DATA_RESIDS)? gretl_VAR_get_residual_matrix(var) :
		NULL;
	} else {
	    M = (code == SYS_DATA_RESIDS)? sys->E : sys->yhat;
	}

	if (M == NULL) {
	    err = E_DATA;
	} else {
	    k = gretl_matrix_cols(M);
	    heads = malloc(k * sizeof *heads);
	    if (heads == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err) {
	    int i, v;

	    for (i=0; i<k && !err; i++) {
		v = (var != NULL)? gretl_VAR_get_variable_number(var, i) :
		    sys->lists[i][1];
		if (v < 0 || v >= dataset->v) {
		    err = E_DATA;
		} else {
		    heads[i] = dataset->varname[v];
		}
	    }
	}

	if (!err) {
	    title = (code == SYS_DATA_RESIDS)? titles[0] : titles[1];
	    wtitle = g_strdup_printf("gretl: %s", _(title));
	    gretl_matrix_print_with_col_heads(M, _(title), heads,
					      dataset, prn);
	}

	free(heads);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	/* FIXME: add matrix as saveable data */
	view_buffer(prn, 80, 400, wtitle, PRINT, NULL);
    }

    g_free(wtitle);
}

static void add_x12_output_menu_item (windata_t *vwin)
{
    const gchar *mpath = "/menubar/Analysis";
    GtkActionEntry entry;

    vwin_menu_add_separator(vwin, mpath);

    action_entry_init(&entry);
    entry.name = "x12aout";
    entry.label = _("View X-12-ARIMA output");
    entry.callback = G_CALLBACK(x12_output_callback);
    vwin_menu_add_item(vwin, mpath, &entry);
}

#include "up_down.h" /* arrows for buttons below */

static GtkWidget *up_down_button (int up)
{
    GtkWidget *img, *w = gtk_button_new();
    GdkPixbuf *pbuf;

    if (up) {
	pbuf = gdk_pixbuf_new_from_inline(-1, up_pixbuf, FALSE, NULL);
    } else {
	pbuf = gdk_pixbuf_new_from_inline(-1, down_pixbuf, FALSE, NULL);
    }

    img = gtk_image_new_from_pixbuf(pbuf);
    gtk_container_add(GTK_CONTAINER(w), img);
    g_object_unref(pbuf);

    return w;
}

static void set_order_vec (GtkWidget *view)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gretl_matrix *m;
    int v, i = 0;

    m = g_object_get_data(G_OBJECT(view), "ordvec");
    model = gtk_tree_view_get_model(GTK_TREE_VIEW(view));

    gtk_tree_model_get_iter_first(model, &iter);
    gtk_tree_model_get(model, &iter, 1, &v, -1);
    m->val[0] = v;

    while (gtk_tree_model_iter_next(model, &iter)) {
	gtk_tree_model_get(model, &iter, 1, &v, -1);
	m->val[++i] = v;
    }
}

static void sensitize_up_down (GtkTreeSelection *selection,
			       GtkWidget *view)
{
    gretl_matrix *v;
    GtkTreeModel *model;
    GtkTreeIter iter;
    GtkTreePath *path;
    GtkWidget *up, *down;
    gboolean s;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(view));
    v = g_object_get_data(G_OBJECT(view), "ordvec");
    up = g_object_get_data(G_OBJECT(view), "up-button");
    down = g_object_get_data(G_OBJECT(view), "down-button");

    gtk_tree_model_get_iter_first(model, &iter);
    s = gtk_tree_selection_iter_is_selected(selection, &iter);
    gtk_widget_set_sensitive(up, !s);

    path = gtk_tree_path_new_from_indices(gretl_vector_get_length(v) - 1,
					  -1);
    s = gtk_tree_selection_path_is_selected(selection, path);
    gtk_widget_set_sensitive(down, !s);
    gtk_tree_path_free(path);
}

static void shift_var_up (GtkButton *b, GtkWidget *view)
{
    GtkTreeSelection *selection;
    GtkTreeModel *model;
    GtkTreeIter seliter, previter;
    GtkTreePath *path;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_get_selected(selection, &model, &seliter);
    path = gtk_tree_model_get_path(model, &seliter);
    if (gtk_tree_path_prev(path)) {
	gtk_tree_model_get_iter(model, &previter, path);
	gtk_list_store_swap(GTK_LIST_STORE(model), &seliter, &previter);
	set_order_vec(view);
	sensitize_up_down(selection, view);
    }
    gtk_tree_path_free(path);
}

static void shift_var_down (GtkButton *b, GtkWidget *view)
{
    GtkTreeSelection *selection;
    GtkTreeModel *model;
    GtkTreeIter seliter, nextiter;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_get_selected(selection, &model, &seliter);
    nextiter = seliter;
    if (gtk_tree_model_iter_next(model, &nextiter)) {
	gtk_list_store_swap(GTK_LIST_STORE(model), &seliter, &nextiter);
	set_order_vec(view);
	sensitize_up_down(selection, view);
    }
}

static void dialog_add_order_selector (GtkWidget *dlg, GRETL_VAR *var,
				       gretl_matrix *ordvec)
{
    GtkWidget *b1, *b2;
    GtkWidget *vbox, *hbox;
    GtkWidget *bbox, *scroller;
    GtkTreeViewColumn *column;
    GtkCellRenderer *renderer;
    GtkTreeSelection *select;
    GtkListStore *store;
    GtkWidget *view, *lbl;
    GtkTreeIter iter;
    const char *vname;
    int i, j, v;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    hbox = gtk_hbox_new(FALSE, 5);
    lbl = gtk_label_new(_("Cholesky ordering:"));
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    store = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_set_data(G_OBJECT(view), "ordvec", ordvec);

    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes(NULL,
						      renderer,
						      "text", 0,
						      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);

    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);
    gtk_tree_view_set_reorderable(GTK_TREE_VIEW(view), FALSE);

    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<var->neqns; i++) {
	j = gretl_vector_get(ordvec, i);
	v = var->ylist[j+1];
	vname = dataset->varname[v];
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, vname,
			   1, j, -1);
    }

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    gtk_tree_selection_select_iter(select, &iter);

    gtk_widget_set_size_request(view, 140 * gui_scale, -1);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW (scroller),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(scroller), view);

#if GTK_MAJOR_VERSION >= 3
    gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(scroller),
					      140 * gui_scale);
#endif

    bbox = gtk_vbox_new(FALSE, 5);
    b1 = up_down_button(1);
    b2 = up_down_button(0);
    gtk_box_pack_start(GTK_BOX(bbox), b1, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(bbox), b2, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(b1), "clicked",
		     G_CALLBACK(shift_var_up), view);
    g_signal_connect(G_OBJECT(b2), "clicked",
		     G_CALLBACK(shift_var_down), view);

    /* FIXME this is ignored in GTK3 ? */
    gtk_widget_set_sensitive(b1, FALSE);
    g_object_set_data(G_OBJECT(view), "up-button", b1);
    g_object_set_data(G_OBJECT(view), "down-button", b2);
    g_signal_connect(G_OBJECT(select), "changed",
		     G_CALLBACK(sensitize_up_down), view);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), scroller, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), bbox, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    g_object_unref(store);
}

static int
impulse_response_setup (GRETL_VAR *var, gretl_matrix *ordvec, int *horizon,
			int *bootstrap, double *alpha, int *piters,
			int *error_bars, GtkWidget *parent)
{
    gchar *title;
    int h = default_VAR_horizon(dataset);
    const char *impulse_opts[] = {
	N_("include bootstrap confidence interval")
    };
    static int active[] = {0};
    GtkWidget *dlg;
    double conf = 1 - *alpha;
    int iters = *piters;
    int resp = -1;

    title = g_strdup_printf("gretl: %s", _("impulse responses"));
    dlg = build_checks_dialog(title, NULL,
			      impulse_opts, 1, active, 0, 0, /* check */
			      0, NULL, /* no radios */
			      &h, _("forecast horizon (periods):"),
			      2, dataset->n / 2, IRF_BOOT,
			      parent, &resp);
    g_free(title);

    if (dlg == NULL) {
	return -1;
    }

    /* includes OPT_E to specify error bars */
    dialog_add_confidence_selector(dlg, &conf, error_bars);

    if (iters == 0) {
	iters = libset_get_int(BOOT_ITERS);
    }
    dialog_add_iters_spinner(dlg, &iters);

    if (ordvec != NULL) {
	dialog_add_order_selector(dlg, var, ordvec);
    }
    gtk_widget_show_all(dlg);

    if (resp < 0) {
	/* canceled */
	*horizon = 0;
    } else {
	*horizon = h;
	*bootstrap = (active[0] > 0);
	if (*bootstrap) {
	    *alpha = 1 - conf;
	    *piters = iters;
	}
    }

    return resp;
}

static int FEVD_setup (GRETL_VAR *var, gretl_matrix *ordvec,
		       int *horizon, int *histogram,
		       GtkWidget *parent)
{
    const char *opts[] = {
	N_("stacked bar graph"),
	N_("line graph"),
    };
    gchar *title;
    int h = default_VAR_horizon(dataset);
    GtkWidget *dlg;
    int optval = 0;
    int resp = -1;

    title = g_strdup_printf("gretl: %s", _("Forecast variance decomposition"));

    dlg = build_checks_dialog(title, NULL, opts,
			      0, NULL, 0, 0, /* no checks */
			      2, &optval, /* two radio buttons */
			      &h, _("forecast horizon (periods):"),
			      2, dataset->n / 2, 0,
			      parent, &resp);
    g_free(title);

    if (dlg == NULL) {
	return -1;
    }

    if (ordvec != NULL) {
	dialog_add_order_selector(dlg, var, ordvec);
    }
    gtk_widget_show_all(dlg);

    if (resp < 0) {
	/* cancelled */
	*horizon = 0;
    } else {
	*horizon = h;
	*histogram = (optval == 0);
    }

    return resp;
}

static void impulse_params_from_action (GtkAction *action,
					int *targ,
					int *shock)
{
    const gchar *s = gtk_action_get_name(action);

    sscanf(s, "Imp:%d:%d", targ, shock);
}

static void FEVD_param_from_action (GtkAction *action,
				    int *targ)
{
    const gchar *s = gtk_action_get_name(action);

    sscanf(s, "FEVD:%d", targ);
}

static int ordvec_default (gretl_matrix *v)
{
    int i, n = gretl_vector_get_length(v);

    for (i=0; i<n; i++) {
	if (v->val[i] != i) {
	    return 0;
	}
    }

    return 1;
}

static gretl_matrix *cholesky_order_vector (GRETL_VAR *var)
{
    gretl_matrix *v = NULL;

    if (var->ord != NULL) {
	v = gretl_matrix_copy(var->ord);
    } else if (var->neqns > 1) {
	int i;

	v = gretl_vector_alloc(var->neqns);
	if (v != NULL) {
	    for (i=0; i<var->neqns; i++) {
		v->val[i] = i;
	    }
	}
    }

    return v;
}

static void FEVD_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int targ, horizon;
    gretl_matrix *ordvec = NULL;
    int histogram = 0;
    int resp, err;

    FEVD_param_from_action(action, &targ);
    ordvec = cholesky_order_vector(var);

    resp = FEVD_setup(var, ordvec, &horizon, &histogram, vwin->main);

    if (resp < 0) {
	/* canceled */
	gretl_matrix_free(ordvec);
	return;
    }

    if (ordvec != NULL) {
	if (ordvec_default(ordvec)) {
	    gretl_matrix_free(ordvec);
	    gretl_VAR_set_ordering(var, NULL);
	} else {
	    gretl_VAR_set_ordering(var, ordvec);
	}
    }

    err = gretl_VAR_plot_FEVD(var, targ, horizon, dataset, histogram);
    gui_graph_handler(err);
}

static void impulse_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int horizon, bootstrap;
    gint shock, targ;
    static double alpha = 0.10;
    static int error_bars = 0;
    static int iters = 0;
    gretl_matrix *ordvec = NULL;
    double this_alpha = 0;
    int save_iters;
    int resp, err;

    impulse_params_from_action(action, &targ, &shock);
    ordvec = cholesky_order_vector(var);
    save_iters = libset_get_int(BOOT_ITERS);

    resp = impulse_response_setup(var, ordvec, &horizon, &bootstrap,
				  &alpha, &iters, &error_bars,
				  vwin->main);

    if (resp < 0) {
	/* canceled */
	gretl_matrix_free(ordvec);
	return;
    }

    if (bootstrap) {
	this_alpha = alpha;
    }

    if (ordvec != NULL) {
	if (ordvec_default(ordvec)) {
	    gretl_matrix_free(ordvec);
	    gretl_VAR_set_ordering(var, NULL);
	} else {
	    gretl_VAR_set_ordering(var, ordvec);
	}
    }

    if (iters != save_iters) {
	libset_set_int(BOOT_ITERS, iters);
    }
    err = gretl_VAR_plot_impulse_response(var, targ, shock,
					  horizon, this_alpha,
					  dataset, error_bars);
    if (iters != save_iters) {
	libset_set_int(BOOT_ITERS, save_iters);
    }

    gui_graph_handler(err);
}

static void multiple_irf_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int horizon, bootstrap;
    static double alpha = 0.10;
    static int error_bars = 0;
    static int iters = 0;
    gretl_matrix *ordvec = NULL;
    double this_alpha = 0;
    int save_iters;
    int resp, err;

    ordvec = cholesky_order_vector(var);
    save_iters = libset_get_int(BOOT_ITERS);

    resp = impulse_response_setup(var, ordvec, &horizon, &bootstrap,
				  &alpha, &iters, &error_bars,
				  vwin->main);

    if (resp < 0) {
	/* canceled */
	gretl_matrix_free(ordvec);
	return;
    }

    if (bootstrap) {
	this_alpha = alpha;
    }

    if (ordvec != NULL) {
	if (ordvec_default(ordvec)) {
	    gretl_matrix_free(ordvec);
	    gretl_VAR_set_ordering(var, NULL);
	} else {
	    gretl_VAR_set_ordering(var, ordvec);
	}
    }

    if (iters != save_iters) {
	libset_set_int(BOOT_ITERS, iters);
    }
    err = gretl_VAR_plot_multiple_irf(var, horizon, this_alpha,
				      dataset, error_bars);
    if (iters != save_iters) {
	libset_set_int(BOOT_ITERS, save_iters);
    }

    gui_graph_handler(err);
}

static int VAR_model_data_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "VarIrf")) {
	return VAR_IRF;
    } else if (!strcmp(s, "VarDecomp")) {
	return VAR_DECOMP;
    } else {
	return VAR_IRF;
    }
}

static void VAR_model_data_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = vwin->data;
    gretl_matrix *ordvec;
    GtkWidget *dlg;
    gchar *title;
    PRN *prn = NULL;
    int code, h = 0;
    int resp = -1;
    int err;

    if (var == NULL) {
	return;
    }

    code = VAR_model_data_code(action);
    h = default_VAR_horizon(dataset);
    ordvec = cholesky_order_vector(var);

    title = g_strdup_printf("gretl: %s",
			    (code == VAR_IRF)? _("impulse responses") :
			    _("variance decompositions"));

    dlg = build_checks_dialog(title, NULL,
			      NULL, 0, NULL, 0, 0, /* no check-buttons */
			      0, NULL,             /* no radios */
			      &h, _("forecast horizon (periods):"),
			      2, dataset->n / 2,
			      0, vwin->main, &resp);

    if (dlg == NULL) {
	goto bailout;
    }

    if (ordvec != NULL) {
	dialog_add_order_selector(dlg, var, ordvec);
    }

    /* blocks till response is selected */
    gtk_widget_show_all(dlg);

    if (resp < 0 || bufopen(&prn)) {
	/* cancel or fail */
	goto bailout;
    }

    if (ordvec != NULL) {
	if (ordvec_default(ordvec)) {
	    gretl_matrix_free(ordvec);
	    ordvec = NULL;
	}
	gretl_VAR_set_ordering(var, ordvec);
	ordvec = NULL;
    }

    if (code == VAR_IRF) {
	err = gretl_VAR_print_all_impulse_responses(var, dataset, h, prn);
    } else if (code == VAR_DECOMP) {
	err = gretl_VAR_print_all_fcast_decomps(var, dataset, h, prn);
    } else {
	err = 1;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	windata_t *viewer;

	viewer = view_buffer_with_parent(vwin, prn, 80, 400, title,
					 code, NULL);
	/* for use when printing in other formats */
	viewer->active_var = h;
    }

 bailout:

    g_free(title);
    gretl_matrix_free(ordvec);
}

static void system_forecast_callback (GtkAction *action, gpointer p)
{
    static gretlopt gopt = OPT_P;
    windata_t *vwin = (windata_t *) p;
    int ci = vwin->role;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    FITRESID *fr;
    int t1, t2, t2est, yno, resp;
    int premax, pre_n, dyn_ok;
    int static_model = 0;
    gretlopt opt = OPT_NONE;
    double conf = 0.95;
    int i, err = 0;

    varnum_from_action(action, &i);

    if (ci == VAR || ci == VECM) {
	var = (GRETL_VAR *) vwin->data;
	t2est = gretl_VAR_get_t2(var);
	yno = var->ylist[i+1];
    } else if (ci == SYSTEM) {
	sys = (equation_system *) vwin->data;
	t2est = sys->t2;
	yno = sys->ylist[i+1];
	static_model = (sys->order == 0);
    } else {
	return;
    }

    t2 = dataset->n - 1;

    /* if no out-of-sample obs are available, alert the user */
    if (t2 == t2est) {
	err = out_of_sample_info(1, &t2);
	if (err) {
	    return;
	}
	t2 = dataset->n - 1;
    }

    /* max number of pre-forecast obs in "best case" */
    premax = dataset->n - 1;

    /* if there are spare obs available, default to an
       out-of-sample forecast */
    if (t2 > t2est) {
	t1 = t2est + 1;
	pre_n = t2est / 2;
	if (pre_n > 100) {
	    pre_n = 100;
	}
	dyn_ok = !static_model;
    } else {
	if (var != NULL) {
	    t1 = levels_order(var);
	} else {
	    t1 = sys->order;
	}
	pre_n = 0;
	dyn_ok = 0;
    }

    resp = forecast_dialog(t1, t1, &t1,
			   t1, t2, &t2, NULL,
			   0, premax, &pre_n,
			   dyn_ok, &gopt, &conf,
			   NULL, vwin->main);
    if (resp < 0) {
	return;
    }

    if (resp == 1) {
	opt = OPT_D;
    } else if (resp == 2) {
	opt = OPT_S;
    }

    fr = get_system_forecast(vwin->data, ci, i, t1, t2, pre_n,
			     dataset, opt, &err);

    if (err) {
	gui_errmsg(err);
    } else {
	char obs1[OBSLEN];
	char obs2[OBSLEN];
	PRN *prn;

	ntolabel(obs1, t1, dataset);
	ntolabel(obs2, t2, dataset);
	lib_command_sprintf("fcast %s %s %s%s", obs1, obs2, dataset->varname[yno],
			    print_flags(opt, FCAST));
	record_command_verbatim();

	if (bufopen(&prn)) {
	    return;
	}

	fr->alpha = 1 - conf;
	err = text_print_forecast(fr, dataset, gopt, prn);
	gui_graph_handler(err);

	view_buffer(prn, (fr->sderr == NULL)? 50 : 78, 400,
		    _("gretl: forecasts"), FCAST, fr);
    }
}

enum {
    SYS_AUTOCORR_TEST,
    SYS_ARCH_TEST,
    SYS_NORMALITY_TEST,
    SYS_RESTRICT
};

static int sys_test_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "autocorr")) {
	return SYS_AUTOCORR_TEST;
    } else if (!strcmp(s, "ARCH")) {
	return SYS_ARCH_TEST;
    } else if (!strcmp(s, "normtest")) {
	return SYS_NORMALITY_TEST;
    } else if (!strcmp(s, "restrict")) {
	return SYS_RESTRICT;
    } else {
	return SYS_NORMALITY_TEST;
    }
}

static void system_test_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    gchar *title = NULL;
    gchar *cstr = NULL;
    PRN *prn;
    int code, order = 0;
    int err = 0;

    if (bufopen(&prn)) {
	return;
    }

    code = sys_test_code(action);

    if (vwin->role == SYSTEM) {
	sys = (equation_system *) vwin->data;
    } else {
	var = (GRETL_VAR *) vwin->data;
    }

    if (code == SYS_AUTOCORR_TEST || code == SYS_ARCH_TEST) {
	int resp;

	order = default_lag_order(dataset);
	resp = spin_dialog((code == SYS_AUTOCORR_TEST)?
			   _("gretl: autocorrelation") :
			   _("gretl: ARCH test"), NULL,
			   &order, _("Lag order for test:"),
			   1, dataset->n / 2, 0, vwin->main);
	if (canceled(resp)) {
	    gretl_print_destroy(prn);
	    return;
	}
    }

    if (code == SYS_AUTOCORR_TEST) {
	title = g_strdup(_("gretl: autocorrelation"));
	cstr = g_strdup_printf("modtest %d --autocorr", order);
	if (var != NULL) {
	    err = gretl_VAR_autocorrelation_test(var, order,
						 dataset,
						 OPT_NONE,
						 prn);
	} else {
	    err = system_autocorrelation_test(sys, order, OPT_NONE, prn);
	}
    } else if (code == SYS_ARCH_TEST) {
	title = g_strdup(_("gretl: ARCH test"));
	cstr = g_strdup_printf("modtest %d --arch", order);
	if (var != NULL) {
	    err = gretl_VAR_arch_test(var, order, dataset, OPT_NONE, prn);
	} else {
	    err = system_arch_test(sys, order, OPT_NONE, prn);
	}
    } else if (code == SYS_NORMALITY_TEST) {
	title = g_strdup_printf("gretl: %s", _("Test for normality of residual"));
	cstr = g_strdup("modtest --normality");
	if (var != NULL) {
	    err = gretl_VAR_normality_test(var, OPT_NONE, prn);
	} else {
	    err = system_normality_test(sys, OPT_NONE, prn);
	}
    } else {
	err = 1;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	add_command_to_stack(cstr, 0);
	view_buffer(prn, 78, 400, title, PRINT, NULL);
    }

    g_free(title);
    g_free(cstr);
}

static void VAR_roots_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int err;

    err = gretl_VAR_roots_plot(var);
    gui_graph_handler(err);
}

static void combined_EC_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int err;

    err = gretl_VECM_combined_EC_plot(var, dataset);
    gui_graph_handler(err);
}

static int sys_ci_from_action (GtkAction *action, int *eqn)
{
    const gchar *s = gtk_action_get_name(action);
    char cmdword[9];

    if (eqn != NULL) {
	sscanf(s, "residplot_%d", eqn);
    }
    sscanf(s, "%*s %8s", cmdword);
    return gretl_command_number(cmdword);
}

static void system_resid_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    int ci, eqn = 0;
    int err;

    ci = sys_ci_from_action(action, &eqn);
    err = gretl_system_residual_plot(vwin->data, ci, eqn, dataset);
    gui_graph_handler(err);
}

static void system_resid_mplot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    int ci = sys_ci_from_action(action, NULL);
    int err;

    err = gretl_system_residual_mplot(vwin->data, ci, dataset);
    gui_graph_handler(err);
}

static const char *VECM_matrix_name (int idx)
{
    if (idx == M_JALPHA) {
	return "alpha";
    } else if (idx == M_JBETA) {
	return "beta";
    } else if (idx == M_JVBETA) {
	return "var(beta)";
    } else {
	/* not registered! */
	return NULL;
    }
}

static int get_VECM_matrix_idx (GtkAction *action)
{
    const gchar *aname = gtk_action_get_name(action);
    int idx;

    sscanf(aname, "matrix %d", &idx);
    return idx;
}

static void VECM_add_matrix (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    gretl_matrix *m = NULL;
    char vname[VNAMELEN];
    const char *mname;
    gchar *blurb = NULL;
    int vecid = gretl_VECM_id(var);
    int idx, resp, show = 1;
    int err = 0;

    idx = get_VECM_matrix_idx(action);
    mname = VECM_matrix_name(idx);
    if (mname == NULL) {
	/* internal error! */
	gui_errmsg(E_DATA);
	return;
    }

    if (!strcmp(mname, "var(beta)")) {
	sprintf(vname, "jvbeta_%d", vecid);
    } else {
	sprintf(vname, "j%s_%d", mname, vecid);
    }
    blurb = g_strdup_printf("%s (%s) from VECM %d\n"
			    "Name (max. %d characters):",
			    mname, gretl_type_get_name(GRETL_TYPE_MATRIX),
			    vecid, VNAMELEN - 1);
    resp = object_name_entry_dialog(vname, GRETL_TYPE_MATRIX, blurb,
				    &show, vwin->main);
    g_free(blurb);
    if (resp < 0) {
	/* canceled */
	return;
    }

    m = gretl_VAR_get_matrix(var, idx, &err);
    if (!err) {
	err = user_var_add_or_replace(vname, GRETL_TYPE_MATRIX, m);
    }
    if (err) {
	gui_errmsg(err);
    } else if (show) {
	view_session();
    }
}

static void add_system_menu_items (windata_t *vwin, int ci)
{
    GtkActionEntry item;
    const gchar *top = "/menubar";
    const gchar *tests = "/menubar/Tests";
    const gchar *save = "/menubar/Save";
    const gchar *graphs = "/menubar/Graphs";
    const gchar *analysis = "/menubar/Analysis";
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    int neqns, nfc, vtarg, vshock;
    char tmp[VNAMELEN2], istr[VNAMELEN];
    char maj[128], min[64];
    const char *cmdword;
    int i, j;

    if (ci == SYSTEM) {
	sys = (equation_system *) vwin->data;
	neqns = sys->neqns;
	nfc = sys->neqns + sys->nidents;
    } else {
	var = (GRETL_VAR *) vwin->data;
	nfc = neqns = gretl_VAR_get_n_equations(var);
    }

    cmdword = gretl_command_word(ci);
    action_entry_init(&item);

    /* FIXME: the following two tests should really be multivariate.
       At present the autocorrelation test for VARs is multivariate
       but the other cases are not.
    */

    if (dataset_is_time_series(dataset)) {
	/* univariate autocorrelation tests */
	item.name = "autocorr";
	item.label = N_("_Autocorrelation");
	item.callback = G_CALLBACK(system_test_call);
	vwin_menu_add_item(vwin, tests, &item);

	/* univariate ARCH tests */
	item.name = "ARCH";
	item.label = N_("A_RCH");
	vwin_menu_add_item(vwin, tests, &item);
    }

    /* multivariate normality test */
    item.name = "normtest";
    item.label = N_("_Normality of residuals");
    item.callback = G_CALLBACK(system_test_call);
    vwin_menu_add_item(vwin, tests, &item);

    if (ci == VECM || ci == SYSTEM) {
	/* linear restrictions */
	item.name = "restrict";
	item.label = N_("Linear restrictions");
	item.callback = G_CALLBACK(gretl_callback);
	vwin_menu_add_item(vwin, tests, &item);
    } else if (ci == VAR) {
	/* regular VAR: omit exogenous variables test */
	if (gretl_VAR_get_exo_list(var) != NULL) {
	    item.name = "VarOmit";
	    item.label = N_("Omit exogenous variables...");
	    item.callback = G_CALLBACK(selector_callback);
	    vwin_menu_add_item(vwin, tests, &item);
	}
	if (var->detflags & DET_TREND) {
	    item.name = "VarOmitTrend";
	    item.label = N_("Omit time trend");
	    item.callback = G_CALLBACK(VAR_omit_auto);
	    vwin_menu_add_item(vwin, tests, &item);
	}
	if (var->detflags & DET_SEAS) {
	    item.name = "VarOmitSeas";
	    item.label = N_("Omit seasonal dummies");
	    item.callback = G_CALLBACK(VAR_omit_auto);
	    vwin_menu_add_item(vwin, tests, &item);
	}
    }

    /* Save residuals */
    for (i=0; i<neqns; i++) {
	sprintf(istr, "resid %d", i);
	sprintf(maj, "%s %d", _("Residuals from equation"), i + 1);
	item.name = istr;
	item.label = maj;
	item.callback = G_CALLBACK(add_system_resid);
	vwin_menu_add_item(vwin, save, &item);
    }

    /* Display residual matrix */
    item.name = "uhat";
    item.label = N_("Display residuals, all equations");
    item.callback = G_CALLBACK(system_data_callback);
    vwin_menu_add_item(vwin, analysis, &item);

    if (ci == SYSTEM) {
	/* Display fitted values matrix */
	item.name = "yhat";
	item.label = N_("Display fitted values, all equations");
	vwin_menu_add_item(vwin, analysis, &item);
    }

    if (neqns > 1) {
	/* Display VCV matrix */
	item.name = "sigma";
	item.label = N_("Cross-equation covariance matrix");
	vwin_menu_add_item(vwin, analysis, &item);
    }

    if (ci == VAR || ci == VECM) {
	/* impulse response printout */
	item.name = "VarIrf";
	item.label = N_("Impulse responses");
	item.callback = G_CALLBACK(VAR_model_data_callback);
	vwin_menu_add_item(vwin, analysis, &item);

	/* variance decomp printout */
	item.name = "VarDecomp";
	item.label = N_("Forecast variance decomposition");
	vwin_menu_add_item(vwin, analysis, &item);
    }

    /* Residual plots */

    action_entry_init(&item);
    item.name = "ResidsMenu";
    item.label = _("Residuals");
    vwin_menu_add_menu(vwin, graphs, &item);

    if (neqns > 1) {
	/* combined residual plot */
	sprintf(min, "comboresid %s", cmdword);
	item.name = min;
	item.label = N_("Combined plot");
	item.callback = G_CALLBACK(system_resid_plot_call);
	vwin_menu_add_item(vwin, "/menubar/Graphs/ResidsMenu", &item);
    }

    if (neqns > 1 && neqns <= 6) {
	/* multiple residual plots in one frame */
	sprintf(min, "multiresid %s", cmdword);
	item.name = min;
	item.label = N_("Multiple plots");
	item.callback = G_CALLBACK(system_resid_mplot_call);
	vwin_menu_add_item(vwin, "/menubar/Graphs/ResidsMenu", &item);
    }

    item.callback = G_CALLBACK(system_resid_plot_call);

    for (i=0; i<neqns; i++) {
	sprintf(min, "residplot_%d %s", i+1, cmdword);
	if (var != NULL) {
	    vtarg = gretl_VAR_get_variable_number(var, i);
	} else {
	    vtarg = sys->ylist[i+1];
	}
	double_underscores(tmp, dataset->varname[vtarg]);
	strcpy(maj, tmp);
	item.name = min;
	item.label = maj;
	vwin_menu_add_item(vwin, "/menubar/Graphs/ResidsMenu", &item);
    }

    /* end residual plots */

    if (ci == VECM) {
	int r = gretl_VECM_rank(var);

	if (r == 1) {
	    item.name = "ecplot";
	    item.label = N_("EC plot");
	    item.callback = G_CALLBACK(combined_EC_plot_call);
	    vwin_menu_add_item(vwin, graphs, &item);
	} else {
	    item.name = "ecplot";
	    item.label = N_("Combined EC plot");
	    item.callback = G_CALLBACK(combined_EC_plot_call);
	    vwin_menu_add_item(vwin, graphs, &item);
	}
    }

    if (ci != SYSTEM) {
	/* VAR inverse roots */
	item.name = "VarRoots";
	item.label = N_("VAR inverse roots");
	item.callback = G_CALLBACK(VAR_roots_plot_call);
	vwin_menu_add_item(vwin, graphs, &item);
    }

    if (ci != SYSTEM && neqns > 1 && neqns <= 4) {
	/* Multiple IRFs */
	item.name = "MultiIrf";
	item.label = N_("Impulse responses (combined)");
	item.callback = G_CALLBACK(multiple_irf_plot_call);
	vwin_menu_add_item(vwin, graphs, &item);
    }

    for (i=0; i<nfc; i++) {
	char newpath[64];
	int dv;

	/* forecast items */
	if (var != NULL) {
	    dv = gretl_VAR_get_variable_number(var, i);
	} else {
	    dv = sys->ylist[i+1];
	}
	double_underscores(tmp, dataset->varname[dv]);
	sprintf(istr, "fcast %d", i);
	item.name = istr;
	item.label = tmp;
	item.callback = G_CALLBACK(system_forecast_callback);
	vwin_menu_add_item(vwin, "/menubar/Analysis/Forecasts", &item);

	if (var == NULL) {
	    continue;
	}

	/* impulse response plots: make menu for target */
	vtarg = gretl_VAR_get_variable_number(var, i);
	double_underscores(tmp, dataset->varname[vtarg]);
	sprintf(istr, "targ_%d", i);
	sprintf(maj, _("Response of %s"), tmp);
	item.name = istr;
	item.label = maj;
	item.callback = NULL;
	vwin_menu_add_menu(vwin, graphs, &item);

	/* path under which to add shocks */
	sprintf(newpath, "/menubar/Graphs/targ_%d", i);

	for (j=0; j<neqns; j++) {
	    /* impulse responses: subitems for shocks */
	    vshock = gretl_VAR_get_variable_number(var, j);
	    double_underscores(tmp, dataset->varname[vshock]);
	    sprintf(istr, "Imp:%d:%d", i, j);
	    sprintf(min, _("to %s"), tmp);
	    item.name = istr;
	    item.label = min;
	    item.callback = G_CALLBACK(impulse_plot_call);
	    vwin_menu_add_item(vwin, newpath, &item);
	}
    }

    if (var != NULL) {
	item.name = "FEVD";
	item.label = _("Forecast variance decomposition");
	item.callback = NULL;
	vwin_menu_add_menu(vwin, graphs, &item);
	for (j=0; j<neqns; j++) {
	    /* FEVD graphs per equation */
	    vtarg = gretl_VAR_get_variable_number(var, j);
	    double_underscores(tmp, dataset->varname[vtarg]);
	    sprintf(istr, "FEVD:%d", j);
	    item.name = istr;
	    item.label = tmp;
	    item.callback = G_CALLBACK(FEVD_plot_call);
	    vwin_menu_add_item(vwin, "/menubar/Graphs/FEVD", &item);
	}
    }

    if (ci == VECM) {
	/* saving things specific to VECMs */
	int mtypes[] = {M_JALPHA, M_JBETA, M_JVBETA};

	/* save error correction terms as series */
	for (i=0; i<jrank(var); i++) {
	    sprintf(istr, "EC %d", i);
	    sprintf(maj, "%s %d", _("EC term"), i+1);
	    item.name = istr;
	    item.label = maj;
	    item.callback = G_CALLBACK(VECM_add_EC_data);
	    vwin_menu_add_item(vwin, save, &item);
	}
	/* save relevant matrices */
	for (i=0; i<G_N_ELEMENTS(mtypes); i++) {
	    sprintf(istr, "matrix %d", mtypes[i]);
	    item.name = istr;
	    item.label = VECM_matrix_name(mtypes[i]);
	    item.callback = G_CALLBACK(VECM_add_matrix);
	    vwin_menu_add_item(vwin, save, &item);
	}
    }

    maybe_add_packages_to_model_menus(vwin);

    if (latex_is_ok()) {
	int n = G_N_ELEMENTS(sys_tex_items);

	vwin_menu_add_menu(vwin, top, &sys_tex_items[0]);
	vwin_menu_add_items(vwin, "/menubar/LaTeX",
			    sys_tex_items + 1, n - 1);
    }
}

void add_system_ui_to_vwin (windata_t *vwin)
{
    vwin_add_ui(vwin, system_items, n_system_items, sys_ui);
    set_model_save_state(vwin, !is_session_model(vwin->data));
    add_system_menu_items(vwin, vwin->role);
    vwin_pack_toolbar(vwin);
    if (vwin->role == VAR || vwin->role == VECM) {
	g_signal_connect(G_OBJECT(vwin->mbar), "button-press-event",
			 G_CALLBACK(check_VAR_menu), vwin);
    }
    gretl_object_ref(vwin->data, (vwin->role == SYSTEM)?
		     GRETL_OBJ_SYS : GRETL_OBJ_VAR);
}

GtkWidget *make_bundle_save_menu (windata_t *vwin)
{
    gretl_bundle *bundle = vwin->data;
    GtkWidget *menu = gtk_menu_new();
    GtkAction *action;
    GtkWidget *item;

    action = gtk_action_new("SaveAs", _("_Save text..."),
			    NULL, NULL);
    g_signal_connect(G_OBJECT(action), "activate",
		     G_CALLBACK(model_output_save), vwin);
    item = gtk_action_create_menu_item(action);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);

    action = gtk_action_new("SaveAsIcon", _("Save bundle to session as _icon"),
			    NULL, NULL);
    g_signal_connect(G_OBJECT(action), "activate",
		     G_CALLBACK(bundle_add_as_icon), vwin);
    item = gtk_action_create_menu_item(action);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    vwin_record_action(vwin, action);

    if (get_user_var_by_data(bundle)) {
	gtk_action_set_sensitive(action, FALSE);
    }

    return menu;
}

static int set_sample_from_model (void *ptr, int role)
{
    MODEL *pmod = NULL;
    int range_set = 0;
    int err = 0;

    if (role == VIEW_MODEL) {
	/* called from single-eqn model window */
	pmod = ptr;
    }

    /* first restore the full dataset */
    err = restore_full_sample(dataset, NULL);

    /* then, if the model was subsampled, restore the subsample */
    if (!err) {
	if (pmod != NULL && pmod->submask != NULL) {
	    err = restrict_sample_from_mask(pmod->submask, dataset,
					    OPT_NONE);
	    range_set = 1;
	} else {
	    /* VAR, VECM or something */
	    int t1 = 0, t2 = 0;

	    model_get_t1_t2(ptr, role, &t1, &t2);
	    if (t1 == 0 && t2 == 0) {
		err = E_DATA;
	    } else if (t1 != dataset->t1 || t2 != dataset->t2) {
		dataset->t1 = t1;
		dataset->t2 = t2;
		range_set = 1;
	    }
	}
    }

    if (err) {
	gui_errmsg(err);
    } else {
	if (range_set) {
	    char comment[64];

	    if (pmod != NULL) {
		sprintf(comment, "# restored sample from model %d\n", pmod->ID);
	    } else {
		strcpy(comment, "# restored sample from model\n");
	    }
	    add_command_to_stack(comment, 0);
	} else {
	    lib_command_strcpy("smpl --full");
	    record_command_verbatim();
	}

	sample_related_menu_state();
	mark_session_changed();
	set_sample_label(dataset);
    }

    return err;
}

/* maybe_set_sample_from_model: return TRUE if the problem situation
   (sample mismatch) is successfully handled, else FALSE.
*/

static gboolean maybe_set_sample_from_model (windata_t *vwin)
{
    const char *msg = N_("The model sample differs from the dataset sample,\n"
			 "so some menu options will be disabled.\n\n"
			 "Do you want to restore the sample on which\n"
			 "this model was estimated?");
    int resp, err = 0;

    resp = yes_no_dialog(NULL, _(msg), vwin_toplevel(vwin));

    if (resp == GRETL_NO) {
	return FALSE;
    }

    err = set_sample_from_model(vwin->data, vwin->role);

    return (err == 0);
}

static gint check_model_menu (GtkWidget *w, GdkEventButton *eb,
			      gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = vwin->data;
    GtkAction *action;
    int resampled = 0;
    int between = 0;
    gboolean s, ok = TRUE;

    if (RQ_SPECIAL_MODEL(pmod)) {
	return FALSE;
    }

    if (dataset == NULL || dataset->Z == NULL) {
	flip(vwin->ui, "/menubar/File/SaveAsIcon", FALSE);
	flip(vwin->ui, "/menubar/File/SaveAndClose", FALSE);
	flip(vwin->ui, "/menubar/Edit/Copy", FALSE);
	flip(vwin->ui, "/menubar/Tests", FALSE);
	flip(vwin->ui, "/menubar/Graphs", FALSE);
	flip(vwin->ui, "/menubar/Analysis", FALSE);
	return FALSE;
    }

    if (pmod->ci == MLE || pmod->ci == GMM ||
	pmod->ci == MPOLS || pmod->ci == BIPROBIT) {
	return FALSE;
    }

    if (model_sample_problem(pmod, dataset)) {
	resampled = (pmod->submask == RESAMPLED);
	between = gretl_is_between_model(pmod);
	ok = FALSE;
    }

    /* check current menu state */
    action = gtk_ui_manager_get_action(vwin->ui, "/menubar/Save/uhat");
    s = gtk_action_is_sensitive(action);

    if (s == ok) {
	/* no need to flip state */
	return FALSE;
    }

    if (s && !ok) {
	if (resampled) {
	    infobox(_("The model sample differs from the dataset sample,\n"
		      "so some menu options will be disabled."));
	} else if (!between) {
	    /* give option to restore model sample */
	    ok = maybe_set_sample_from_model(vwin);
	    if (ok) {
		return FALSE;
	    }
	}
    }

    flip(vwin->ui, "/menubar/Analysis/Forecasts", ok);
    flip(vwin->ui, "/menubar/Save/yhat", ok);
    flip(vwin->ui, "/menubar/Save/uhat", ok);
    flip(vwin->ui, "/menubar/Save/uhat2", ok);
    flip(vwin->ui, "/menubar/Save/NewVar", ok);

    if (resampled) {
	flip(vwin->ui, "/menubar/Tests", FALSE);
	flip(vwin->ui, "/menubar/Analysis/DisplayAFR", FALSE);
	flip(vwin->ui, "/menubar/Analysis/Bootstrap", FALSE);
	flip(vwin->ui, "/menubar/Graphs", FALSE);
    }

    if (between) {
	flip(vwin->ui, "/menubar/Edit/Revise", FALSE);
    }

    return FALSE;
}

static gint check_VAR_menu (GtkWidget *w, GdkEventButton *eb,
			    gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    GtkAction *action;
    gboolean s, ok = TRUE;

    if (complex_subsampled()) {
	ok = FALSE;
    }

    action = gtk_ui_manager_get_action(vwin->ui, "/menubar/Tests");
    s = gtk_action_is_sensitive(action);

    if (s == ok) {
	/* no need to flip state */
	return FALSE;
    }

    flip(vwin->ui, "/menubar/Edit/Revise", ok);
    flip(vwin->ui, "/menubar/Tests", ok);
    flip(vwin->ui, "/menubar/Save", ok);
    flip(vwin->ui, "/menubar/Graphs", ok);
    flip(vwin->ui, "/menubar/Analysis/Forecasts", ok);
    flip(vwin->ui, "/menubar/Analysis/VarIrf", ok);
    flip(vwin->ui, "/menubar/Analysis/VarDecomp", ok);

    if (!ok) {
	warnbox(_("dataset is subsampled"));
    }

    return FALSE;
}

static int object_overwrite_ok (const char *name, GretlType t,
				GtkWidget *parent)
{
    gchar *info = name_conflict_message(name, t);
    gchar *msg = g_strdup_printf("%s\n%s", info, _("OK to overwrite it?"));
    int resp;

    resp = yes_no_dialog("gretl", msg, parent);
    g_free(info);
    g_free(msg);

    return (resp == GRETL_YES);
}

/* note: returns non-zero if the varname is not acceptable */

static int real_gui_validate_varname (const char *name,
				      GretlType t,
				      int allow_overwrite,
				      GtkWidget *parent)
{
    int i, n = strlen(name);
    char namebit[VNAMELEN];
    unsigned char c;
    int err = 0;

    *namebit = '\0';

    if (n > VNAMELEN - 1) {
	strncat(namebit, name, VNAMELEN - 1);
	errbox_printf(_("Variable name %s... is too long\n"
			"(the max is %d characters)"), namebit,
		      VNAMELEN - 1);
	err = 1;
    } else if (!(isalpha(*name))) {
	errbox_printf(_("First char of name ('%c') is bad\n"
			"(first must be alphabetical)"), *name);
	err = 1;
    } else {
	for (i=1; i<n && !err; i++) {
	    c = (unsigned char) name[i];

	    if ((!(isalpha(c)) && !(isdigit(c)) && c != '_') || c > 127) {
		errbox_printf(_("Name contains an illegal char (in place %d)\n"
				"Use only unaccented letters, digits and underscore"),
			      i + 1);
		err = 1;
	    }
	}
    }

    if (!err && t != GRETL_TYPE_NONE) {
	/* check for variable type collisions */
	GretlType t0 = gretl_type_from_name(name, dataset);

	if (t0 != GRETL_TYPE_NONE) {
	    /* there's already a variable of this name */
	    if (t == t0 && allow_overwrite) {
		/* the types agree: overwrite? */
		if (allow_overwrite < 2) {
		    err = !object_overwrite_ok(name, t, parent);
		}
	    } else {
		/* the types disgree: won't work */
		gchar *msg = name_conflict_message(name, t0);

		errbox(msg);
		g_free(msg);
		err = 1;
	    }
	}
    }

    return err;
}

/* The "gui_validate_varname" family: the functions below check the
   putative @name for legality as a gretl variable name and return
   non-zero if it's not legal. In addition, both return non-zero if
   the name is valid but belongs to an existing variable of a type
   other than @type.

   They diverge in this respect:

   gui_validate_varname_strict: unconditionally returns non-zero
   if there's an existing variable of the same name, even if it's
   of type @type.

   gui_validate_varname: checks with the user whether overwriting is
   OK in the case where a variable of type @type already exists; if
   so, it's assumed that the distinction between redefining a variable
   and adding a new variable is handled by the caller.

   gui_validate_varname_easy: assumes that overwriting an existing
   variable of suitable type is OK.
*/

int gui_validate_varname_strict (const char *name, GretlType type,
				 GtkWidget *parent)
{
    return real_gui_validate_varname(name, type, 0, parent);
}

int gui_validate_varname (const char *name, GretlType type,
			  GtkWidget *parent)
{
    return real_gui_validate_varname(name, type, 1, parent);
}

int gui_validate_varname_easy (const char *name, GretlType type)
{
    return real_gui_validate_varname(name, type, 2, NULL);
}

gint popup_menu_handler (GtkWidget *widget, GdkEventButton *event,
			 gpointer data)
{
    if (right_click(event)) {
	gtk_menu_popup(GTK_MENU(data), NULL, NULL, NULL, NULL,
		       event->button, event->time);
	return TRUE;
    }

    return FALSE;
}

void add_popup_item (const gchar *label, GtkWidget *menu,
		     GCallback callback,
		     gpointer data)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(label);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    g_signal_connect(G_OBJECT(item), "activate",
		     G_CALLBACK(callback), data);
    gtk_widget_show(item);
}
