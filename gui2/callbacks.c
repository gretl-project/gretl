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

#include "boxplots.h"

/* these live in dialogs.c */
extern GtkWidget *active_edit_id;
extern GtkWidget *active_edit_name;
extern GtkWidget *active_edit_text;

static void doubleclick_action (windata_t *vwin)
{
    switch (vwin->role) {
    case MAINWIN:
	if (datainfo != NULL && datainfo->n > 0) {
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
	db_series_callback(NULL, vwin);
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
    { OPEN_DATA,       "OpenGdt" },
    { APPEND_DATA,     "AppendGdt" },
    { OPEN_ASCII,      "OpenASCII" },
    { APPEND_ASCII,    "AppendASCII" },
    { OPEN_CSV,        "OpenCSV" },
    { APPEND_CSV,      "AppendCSV" },
    { OPEN_OCTAVE,     "OpenOctave" },
    { APPEND_OCTAVE,   "AppendOctave" },
    { OPEN_GNUMERIC,   "OpenGnumeric" },
    { APPEND_GNUMERIC, "AppendGnumeric" },
    { OPEN_XLS,        "OpenXLS" },
    { APPEND_XLS,      "AppendXLS" },
    { OPEN_WF1,        "OpenWF1" },
    { APPEND_WF1,      "AppendWF1" },
    { OPEN_DTA,        "OpenDTA" },
    { APPEND_DTA,      "AppendDTA" },
    { OPEN_SAV,        "OpenSAV" },
    { APPEND_SAV,      "AppendSAV" },
    { OPEN_JMULTI,     "OpenJMulTi" },
    { APPEND_JMULTI,   "AppendJMulTi" },
    { OPEN_ODS,        "OpenODS" },
    { APPEND_ODS,      "AppendODS" },
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
    int code;

    if (dataset_locked()) {
	return;
    }

    code = open_data_code(gtk_action_get_name(action));

    switch (code) {
    case OPEN_DATA:
    case APPEND_DATA:
    case OPEN_ASCII:
    case APPEND_ASCII:
	file_selector(_("Open data file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_CSV:
    case APPEND_CSV:
	if (delimiter_dialog(NULL)) {
	    return;
	}
	file_selector(_("Open CSV file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_OCTAVE:
    case APPEND_OCTAVE:
	file_selector(_("Open Octave file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_GNUMERIC:
    case APPEND_GNUMERIC:
	file_selector(_("Open Gnumeric file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_XLS:
    case APPEND_XLS:
	file_selector(_("Open Excel file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_WF1:
    case APPEND_WF1:
	file_selector(_("Open Eviews workfile"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_DTA:
    case APPEND_DTA:
	file_selector(_("Open Stata file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_SAV:
    case APPEND_SAV:
	file_selector(_("Open SPSS file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_JMULTI:
    case APPEND_JMULTI:
	file_selector(_("Open JMulTi file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_ODS:
    case APPEND_ODS:
	file_selector(_("Open ODS file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_MARKERS:
	file_selector(_("gretl: add markers"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_RATS_DB:
    case OPEN_PCGIVE_DB:
	file_selector(_("gretl: open database"), code, FSEL_DATA_NONE, NULL);
	break;
    default:
	fprintf(stderr, "open_data: unrecognized action '%s'\n",
		gtk_action_get_name(action));
	break;
    }
}

void open_script (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "OpenScript")) {
	file_selector(_("Open script file"), OPEN_SCRIPT, FSEL_DATA_NONE, NULL);
    } else if (!strcmp(s, "OpenSession")) {
	file_selector(_("Open session file"), OPEN_SESSION, FSEL_DATA_NONE, NULL);
    }
}

void file_save (windata_t *vwin, int ci)
{
    gretlopt opt = OPT_NONE;
    gpointer p = NULL;

    switch (ci) {
    case SAVE_OUTPUT:
	file_selector(_("Save output file"), ci, FSEL_DATA_VWIN, vwin);
	break;
    case SAVE_CONSOLE:
	file_selector(_("Save console output"), ci, FSEL_DATA_VWIN, vwin);
	break;
    case SAVE_SCRIPT:
	file_selector(_("Save command script"), ci, FSEL_DATA_VWIN, vwin);
	break;
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case SAVE_DBDATA:
    case SAVE_FUNCTIONS:	
	data_save_selection_wrapper(ci, NULL);
	break;
    case EXPORT_CSV:
	if (delimiter_dialog(&opt)) {
	    return;
	}
	p = GINT_TO_POINTER(opt);
    case EXPORT_R:
    case EXPORT_OCTAVE:
    case EXPORT_DAT:
    case EXPORT_JM:
	data_save_selection_wrapper(ci, p);
	break;
    case SAVE_TEX:
	file_selector(_("Save LaTeX file"), ci, FSEL_DATA_MISC, vwin->data);
	break;
    case SAVE_TEXT:
	file_selector(_("Save text"), ci, FSEL_DATA_MISC, vwin->data);
	break;
    case SAVE_GP_CMDS:
	file_selector(_("Save gnuplot commands"), ci, FSEL_DATA_VWIN, vwin);
	break;
    case SAVE_R_CMDS:
	file_selector(_("Save R commands"), ci, FSEL_DATA_VWIN, vwin);
	break;
    default:
	dummy_call();
    }
}

static int fsave_code (const gchar *s)
{
    if (!strcmp(s, "SaveAsGdt"))
	return SAVE_DATA;
    if (!strcmp(s, "SaveAsDb"))
	return SAVE_DBDATA;
    if (!strcmp(s, "ExportCSV"))
	return EXPORT_CSV;
    if (!strcmp(s, "ExportR"))
	return EXPORT_R;
    if (!strcmp(s, "ExportOctave"))
	return EXPORT_OCTAVE;
    if (!strcmp(s, "ExportPcGive"))
	return EXPORT_DAT;
    if (!strcmp(s, "ExportJMulTi"))
	return EXPORT_JM;
    if (!strcmp(s, "NewGfn"))
	return SAVE_FUNCTIONS;

    return SAVE_DATA;
}

void fsave_callback (GtkAction *action, gpointer p)
{
    const gchar *s = gtk_action_get_name(action);
    int ci = fsave_code(s);

    file_save(p, ci);
}

void dummy_call (void)
{
    errbox(_("Sorry, this item not yet implemented!"));
}

void print_report (GtkAction *action)
{
    PRN *prn;

    if (bufopen(&prn)) return;

    data_report(datainfo, &paths, prn);

    view_buffer(prn, 77, 400, _("gretl: data summary"), 
		DATA_REPORT, NULL);
}

void edit_header (GtkAction *action)
{
    if (data_status & BOOK_DATA) {
	errbox(_("You don't have permission to do this"));
    } else { 
	edit_buffer(&datainfo->descrip, 80, 400, _("gretl: edit data info"),
		    EDIT_HEADER);
    }
}

static int model_action_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    int ci = gretl_command_number(s);

    if (ci == 0) {
	if (!strcmp(s, "CORC"))
	    ci = CORC;
	else if (!strcmp(s, "HILU"))
	    ci = HILU;
	else if (!strcmp(s, "PWE"))
	    ci = PWE;
	else if (!strcmp(s, "PANEL_WLS"))
	    ci = PANEL_WLS;
	else if (!strcmp(s, "PANEL_B"))
	    ci = PANEL_B;
	else if (!strcmp(s, "VLAGSEL"))
	    ci = VLAGSEL;
    }

    return ci;
}

void fit_resid_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *mydata = (windata_t *) data; 
    MODEL *pmod = mydata->data;
    int code = GENR_RESID; 

    if (!strcmp(s, "yhat")) {
	code = GENR_FITTED;
    } else if (!strcmp(s, "uhat")) {
	code = GENR_RESID;
    } else if (!strcmp(s, "uhat2")) {
	code = GENR_RESID2;
    } else if (!strcmp(s, "h")) {
	code = GENR_H;
    } else if (!strcmp(s, "ahat")) {
	code = GENR_AHAT;
    }

    add_fit_resid(pmod, code, 0);
}

void model_stat_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) data; 
    MODEL *pmod = vwin->data;
    int code = ESS; 

    if (!strcmp(s, "ess")) {
	code = ESS;
    } else if (!strcmp(s, "se")) {
	code = SIGMA;
    } else if (!strcmp(s, "rsq")) {
	code = R2;
    } else if (!strcmp(s, "trsq")) {
	code = TR2;
    } else if (!strcmp(s, "df")) {
	code = DF;
    } else if (!strcmp(s, "lnL")) {
	code = LNL;
    } else if (!strcmp(s, "AIC")) {
	code = AIC;
    } else if (!strcmp(s, "BIC")) {
	code = BIC;
    } else if (!strcmp(s, "HQC")) {
	code = HQC;
    }

    add_model_stat(pmod, code);
}

void model_callback (GtkAction *action, gpointer data) 
{
    int code = model_action_code(action);

    modelspec_dialog(code);
}

void model_genr_callback (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    edit_dialog(_("gretl: add var"), 
		_("Enter formula for new variable:"),
		"", do_model_genr, vwin, 
		MODEL_GENR, VARCLICK_INSERT_NAME, NULL);   
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
    if (!strcmp(s, "VLAGSEL"))
	return VLAGSEL;
    if (!strcmp(s, "ConfEllipse"))
	return ELLIPSE;
    if (!strcmp(s, "VarOmit"))
	return VAROMIT;

    return 0;
}

void selector_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) data;
    char title[64];
    int ci;

    ci = selector_callback_code(s);

    if (ci == ADD || ci == OMIT || ci == COEFFSUM || ci == ELLIPSE) {
	set_window_busy(vwin);
    }

    strcpy(title, "gretl: ");

    if (ci == COINT || ci == COINT2) {
	selection_dialog(_("gretl: cointegration test"), do_coint, ci);
    } else if (ci == VAR || ci == VECM) {
	selection_dialog((ci == VAR)? _("gretl: VAR") : _("gretl: VECM"),
			 do_vector_model, ci);
    } else if (ci == VLAGSEL) {
	selection_dialog(_("gretl: VAR lag selection"), do_vector_model, ci);
    } else if (ci == GR_XY || ci == GR_IMP || ci == GR_DUMMY
	       || ci == SCATTERS || ci == GR_3D
	       || ci == GR_XYZ) {
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
	    selfunc = do_scatters;
	    break;
	default:
	    return;
	}
	selection_dialog(_("gretl: define graph"), selfunc, ci);
    } else if (ci == ADD || ci == OMIT) {
	simple_selection(_("gretl: model tests"), do_add_omit, ci, vwin);
    } else if (ci == VAROMIT) {
	simple_selection(_("gretl: model tests"), do_VAR_omit, ci, vwin);
    } else if (ci == COEFFSUM) {
	simple_selection(_("gretl: model tests"), do_coeff_sum, ci, vwin);
    } else if (ci == ELLIPSE) {
	simple_selection(_("gretl: model tests"), do_confidence_region, ci, vwin);
    } else if (ci == GR_PLOT) {
	simple_selection(_("gretl: define graph"), do_graph_from_selector, ci, NULL);
    } else if (ci == TSPLOTS) {
	simple_selection(_("gretl: define graph"), do_scatters, ci, vwin);
    } else if (ci == SPEARMAN) {
	strcat(title, _("rank correlation"));
	simple_selection(title, do_rankcorr, ci, vwin);
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
    if (!strcmp(s, "GR_BOX")) 
	return GR_BOX;
    if (!strcmp(s, "GR_NBOX")) 
	return GR_NBOX;
    if (!strcmp(s, "gmm")) 
	return GMM;
    if (!strcmp(s, "mle")) 
	return MLE;
    if (!strcmp(s, "nls")) 
	return NLS;
    if (!strcmp(s, "system")) 
	return SYSTEM;
    if (!strcmp(s, "restrict")) 
	return RESTRICT;
    if (!strcmp(s, "MINIBUF")) 
	return MINIBUF;
    return 0;
}

void gretl_callback (GtkAction *action, gpointer data)
{
    const char *title = NULL;
    const char *query = NULL;
    const char *defstr = NULL;
    void (*okfunc)() = NULL;
    guint varclick = VARCLICK_NONE;
    int cmd, cancel = 0;

    cmd = gretl_callback_code(gtk_action_get_name(action));

    switch (cmd) {
    case SMPLBOOL:
	title = N_("gretl: restrict sample");
	query = N_("Enter boolean condition for selecting cases:");
	okfunc = do_samplebool;
	varclick = VARCLICK_INSERT_NAME;
	break;
    case GENR:
	title = N_("gretl: add var");
	query = N_("Enter formula for new variable\n"
		   "(or just a name, to enter data manually)");
	okfunc = do_genr;
	varclick = VARCLICK_INSERT_NAME;
	break;
    case VSETMISS:
	title = N_("gretl: missing code");
	query = N_("Enter value to be read as \"missing\":");
	okfunc = do_variable_setmiss;
	break;
    case GSETMISS:
	title = N_("gretl: missing code");
	query = N_("Enter value to be read as \"missing\":");
	okfunc = do_global_setmiss;
	break;
    case GR_BOX:
    case GR_NBOX:
	title = N_("gretl: boxplots");
	query = N_("Specify variables to plot:");
	okfunc = do_box_graph;
	varclick = VARCLICK_INSERT_NAME;
	defstr = get_last_boxplots_string();
	break;
    case GMM:
	title = N_("gretl: GMM");
	query = N_("GMM: Specify function and orthogonality conditions:");
	okfunc = do_gmm_model;
	varclick = VARCLICK_INSERT_TEXT;
	break;	
    case MLE:
	title = N_("gretl: maximum likelihood");
	query = N_("MLE: Specify function, and derivatives if possible:");
	okfunc = do_mle_model;
	varclick = VARCLICK_INSERT_TEXT;
	break;	
    case NLS:
	title = N_("gretl: nonlinear least squares");
	query = N_("NLS: Specify function, and derivatives if possible:");
	okfunc = do_nls_model;
	varclick = VARCLICK_INSERT_TEXT;
	break;	
    case SYSTEM:
	title = N_("gretl: simultaneous equations system");
	query = N_("Specify simultaneous equations:");
	data = NULL;
	okfunc = do_eqn_system;
	varclick = VARCLICK_INSERT_TEXT;
	break;
    case RESTRICT:
	title = N_("gretl: linear restrictions");
	query = N_("Specify restrictions:");
	okfunc = do_restrict;
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

    edit_dialog(_(title), _(query), defstr, okfunc, data, 
		cmd, varclick, &cancel);   
}

void genr_callback (void)
{
    edit_dialog(_("gretl: add var"), 
		_("Enter formula for new variable\n"
		  "(or just a name, to enter data manually)"),
		NULL, do_genr, NULL, 
		GENR, VARCLICK_INSERT_NAME, NULL);   
}

void minibuf_callback (void)
{
    edit_dialog(_("gretl: command entry"),
		_("Type a command:"),
		NULL, do_minibuf, NULL, 
		MINIBUF, VARCLICK_NONE, NULL);   
}

void newdata_callback (void) 
{
    int resp, n = 50;

    if (dataset_locked()) {
	return;
    }

    resp = spin_dialog(_("gretl: create data set"), NULL, &n, 
		       _("Number of observations:"), 
		       2, 100000, 0);

    if (resp < 0) {
	/* canceled */
	return;
    }

    if (open_nulldata(&Z, datainfo, data_status, n, NULL)) {
	errbox(_("Failed to create empty data set"));
	return;
    }

    new_data_structure_dialog();
}

void xcorrgm_callback (void)
{
    if (mdata_selection_count() == 2) {
	do_xcorrgm(NULL);
    } else {
	char title[64];

	strcpy(title, "gretl: ");
	strcat(title, _("cross-correlogram"));
	simple_selection(title, do_xcorrgm, XCORRGM, NULL);
    }
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
    void *handle;
    int (*run_nist_tests)(const char *, const char *, int);
    gchar *datadir = NULL;
    gchar *fname = NULL;
    
    run_nist_tests = gui_get_plugin_function("run_nist_tests", 
					     &handle);
    if (run_nist_tests == NULL) {
	return;
    }

    datadir = g_strdup_printf("%sdata%s", paths.gretldir, SLASHSTR);
    fname = g_strdup_printf("%snist.out", paths.dotdir);

    (*run_nist_tests)(datadir, fname, nist_verbosity(action));

    close_plugin(handle);

    view_file(fname, 0, 1, 78, 400, VIEW_CODEBOOK);

    g_free(datadir);
    g_free(fname);
}

#if defined (ENABLE_MAILER) && !defined(G_OS_WIN32)

void send_file (char *fullname)
{
    int (*email_file) (const char *, const char *);
    void *handle;

    email_file = gui_get_plugin_function("email_file", &handle);
    if (email_file == NULL) {
        return;
    }
    
    email_file(fullname, paths.dotdir);
    close_plugin(handle);
}

#endif
