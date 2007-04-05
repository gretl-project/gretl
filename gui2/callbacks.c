/*
 *   Copyright (C) Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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
#include "boxplots.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "ssheet.h"
#include "treeutils.h"

/* these live in dialogs.c */
extern GtkWidget *active_edit_id;
extern GtkWidget *active_edit_name;
extern GtkWidget *active_edit_text;

static void doubleclick_action (windata_t *win)
{
    switch (win->role) {
    case MAINWIN:
	display_var();
	break;
    case TEXTBOOK_DATA:
	browser_open_data(NULL, win);
	break;
    case PS_FILES:
	browser_open_ps(NULL, win);
	break;
    case FUNC_FILES:
	browser_call_func(NULL, win);
	break;
    case NATIVE_DB:
	open_db_index(NULL, win); 
	break;
    case REMOTE_DB:
	open_remote_db_index(NULL, win);
	break;
    case NATIVE_SERIES:
    case RATS_SERIES:
    case PCGIVE_SERIES:
    case REMOTE_SERIES:
	gui_get_db_series(win, 0, NULL);
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
        gdk_window_get_pointer (event->window, &x, &y, &state);
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
	if (select == NULL) return FALSE;
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

void open_data (gpointer data, guint code, GtkWidget *widget)
{
    if (dataset_locked()) {
	return;
    }

    switch (code) {
    case OPEN_DATA:
    case APPEND_DATA:
    case OPEN_ASCII:
    case APPEND_ASCII:
	file_selector(_("Open data file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_CSV:
    case APPEND_CSV:
	delimiter_dialog(NULL);
	file_selector(_("Open CSV file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_OCTAVE:
    case APPEND_OCTAVE:
	file_selector(_("Open Octave file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_BOX:
	file_selector(_("Open BOX file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_GNUMERIC:
    case APPEND_GNUMERIC:
	file_selector(_("Open Gnumeric file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_EXCEL:
    case APPEND_EXCEL:
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
    case OPEN_JMULTI:
    case APPEND_JMULTI:
	file_selector(_("Open JMulTi file"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_MARKERS:
	file_selector(_("gretl: add markers"), code, FSEL_DATA_NONE, NULL);
	break;
    case OPEN_RATS_DB:
    case OPEN_PCGIVE_DB:
	file_selector(_("gretl: open database"), code, FSEL_DATA_NONE, NULL);
	break;
    default:
	errbox("Unrecognized data code");
	break;
    }
}

void open_script (gpointer data, guint action, GtkWidget *widget)
{
    if (action == OPEN_SCRIPT) {
	file_selector(_("Open script file"), action, FSEL_DATA_NONE, NULL);
    } else if (action == OPEN_SESSION) {
	file_selector(_("Open session file"), action, FSEL_DATA_NONE, NULL);
    }
}

void file_save (gpointer data, guint file_code, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    gretlopt opt = OPT_NONE;
    gpointer p = NULL;

    switch (file_code) {
    case SAVE_OUTPUT:
	file_selector(_("Save output file"), SAVE_OUTPUT, FSEL_DATA_MISC, vwin);
	break;
    case SAVE_CONSOLE:
	file_selector(_("Save console output"), SAVE_CONSOLE, FSEL_DATA_MISC, vwin);
	break;
    case SAVE_SCRIPT:
	file_selector(_("Save command script"), SAVE_SCRIPT, FSEL_DATA_MISC, vwin);
	break;
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case SAVE_DBDATA:
    case SAVE_FUNCTIONS:	
	data_save_selection_wrapper(file_code, NULL);
	break;
    case EXPORT_CSV:
	delimiter_dialog(&opt);
	p = GINT_TO_POINTER(opt);
    case EXPORT_R:
    case EXPORT_OCTAVE:
    case EXPORT_DAT:
    case EXPORT_JM:
	data_save_selection_wrapper(file_code, p);
	break;
    case SAVE_TEX:
	file_selector(_("Save LaTeX file"), file_code, FSEL_DATA_MISC, vwin->data);
	break;
    case SAVE_GP_CMDS:
	file_selector(_("Save gnuplot commands"), file_code, FSEL_DATA_MISC, vwin);
	break;
    default:
	dummy_call();
    }
}

void dummy_call (void)
{
    errbox(_("Sorry, this item not yet implemented!"));
}

void print_report (gpointer data, guint unused, GtkWidget *widget)
{
    PRN *prn;

    if (bufopen(&prn)) return;

    data_report(datainfo, &paths, prn);

    view_buffer(prn, 77, 400, _("gretl: data summary"), 
		DATA_REPORT, NULL);
}

void edit_header (gpointer data, guint unused, GtkWidget *widget)
{
    if (data_status & BOOK_DATA) {
	errbox(_("You don't have permission to do this"));
    } else { 
	edit_buffer(&datainfo->descrip, 80, 400, _("gretl: edit data info"),
		    EDIT_HEADER);
    }
}

void fit_resid_callback (gpointer data, guint code, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data; 
    MODEL *pmod = mydata->data;

    add_fit_resid(pmod, code, 0);
}

void VAR_resid_callback (gpointer data, guint eqnum, GtkWidget *widget)
{
    add_system_resid(data, eqnum, VAR);
}

void SYS_resid_callback (gpointer data, guint eqnum, GtkWidget *widget)
{
    add_system_resid(data, eqnum, SYSTEM);
}

void model_stat_callback (gpointer data, guint which, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data; 
    MODEL *pmod = mydata->data;

    add_model_stat(pmod, which);
}

void model_callback (gpointer data, guint model_code, GtkWidget *widget) 
{
    int presel = 0;

    if (widget == NULL && data != NULL) {
	/* preselected dependent variable */
	presel = GPOINTER_TO_INT(data);
    }

    selection_dialog(_("gretl: specify model"), do_model, model_code,
		     presel);
}

void model_genr_callback (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data;

    edit_dialog(_("gretl: add var"), _("Enter formula for new variable:"),
		"", do_model_genr, mydata, 
		MODEL_GENR, VARCLICK_INSERT_NAME, NULL);   
}

void selector_callback (gpointer data, guint action, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    char title[64];

    if (action == ADD || action == OMIT || action == COEFFSUM ||
	action == ELLIPSE) {
	set_window_busy(vwin);
    }

    strcpy(title, "gretl: ");

    if (action == COINT || action == COINT2) {
	selection_dialog(_("gretl: cointegration test"), do_coint, action, 0);
    } else if (action == VAR || action == VECM) {
	selection_dialog((action == VAR)? _("gretl: VAR") : _("gretl: VECM"),
			 do_vector_model, action, 0);
    } else if (action == VLAGSEL) {
	selection_dialog(_("gretl: VAR lag selection"), do_vector_model, action, 0);
    } else if (action == GR_XY || action == GR_IMP || action == GR_DUMMY
	       || action == SCATTERS || action == GR_3D) {
	int (*selfunc)() = NULL;

	switch (action) {
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
	case SCATTERS:
	    selfunc = do_scatters;
	    break;
	default:
	    return;
	}
	selection_dialog(_("gretl: define graph"), selfunc, action, 0);
    } else if (action == ADD || action == OMIT) {
	simple_selection(_("gretl: model tests"), do_add_omit, action, vwin);
    } else if (action == VAROMIT) {
	simple_selection(_("gretl: model tests"), do_VAR_omit, action, vwin);
    } else if (action == COEFFSUM) {
	simple_selection(_("gretl: model tests"), do_coeff_sum, action, vwin);
    } else if (action == ELLIPSE) {
	simple_selection(_("gretl: model tests"), do_confidence_region, action, vwin);
    } else if (action == GR_PLOT) {
	simple_selection(_("gretl: define graph"), do_graph_from_selector, action, NULL);
    } else if (action == TSPLOTS) {
	simple_selection(_("gretl: define graph"), do_scatters, action, vwin);
    } else if (action == SPEARMAN) {
	strcat(title, _("rank correlation"));
	simple_selection(title, do_spearman, action, vwin);
    } else {
	errbox("selector_callback: code was not recognized");
    }
}

void gretl_callback (gpointer data, guint action, GtkWidget *widget)
{
    const char *title = NULL;
    const char *query = NULL;
    const char *defstr = NULL;
    void (*okfunc)() = NULL;
    guint varclick = VARCLICK_NONE;

    switch (action) {
    case SMPLBOOL:
	title = N_("gretl: restrict sample");
	query = N_("Enter boolean condition for selecting cases:");
	okfunc = do_samplebool;
	varclick = VARCLICK_INSERT_NAME;
	break;
    case GENR:
	title = N_("gretl: add var");
	query = N_("Enter formula for new variable:");
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
	defstr = get_boxplots_string();
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
	errbox("Bug: unrecognized action code in gretl_callback");
	return;
    }

    edit_dialog(_(title), _(query), defstr, okfunc, data, 
		action, varclick, NULL);   
}

void gp_send_callback (GtkWidget *w, gpointer data)
{
    gp_to_gnuplot(data, 0, w);
}

void file_save_callback (GtkWidget *w, gpointer data)
{
    guint u = 0;
    windata_t *vwin = (windata_t *) data;

    if (g_object_get_data(G_OBJECT(vwin->dialog), "text_out")) {
	u = SAVE_OUTPUT;
    } else {
	switch (vwin->role) {
	case EDIT_SCRIPT:
	case VIEW_SCRIPT:
	case VIEW_LOG:
	case VIEW_FUNC_CODE:
	    u = SAVE_SCRIPT;
	    break;
	case GR_PLOT:
	    u = SAVE_GP_CMDS;
	    break;
	default:
	    errbox(_("Sorry, not yet implemented"));
	    return;
	}
    }

    file_save(data, u, w);
}

void add_rand_callback (gpointer data, guint i, GtkWidget *widget) 
{
    const char *title[] = {
	N_("gretl: uniform variable"),
	N_("gretl: normal variable"),
	N_("gretl: chi-square variable"),
	N_("gretl: Student's t variable"),
	N_("gretl: binomial variable"),
	N_("gretl: Poisson variable")
    };
    const char *info[] = {
	N_("Enter name for variable, and\n"
	   "minimum and maximum values:"), 
	N_("Enter name, mean and standard deviation:"), 
	N_("Enter name and degrees of freedom:"), 
	N_("Enter name and degrees of freedom:"), 
	N_("Enter name, number of trials and probability:"), 
	N_("Enter name and mean (scalar or series):")
    };
    const char *template[] = {
	"unif 0 1",
	"norm 0 1", 
	"chi 5", 
	"st 20", 
	"bin 10 0.5", 
	"pois xbar"
    };
	
    edit_dialog(_(title[i]), _(info[i]), _(template[i]),
		add_rand_series, GUINT_TO_POINTER(i), 
		GENR_RANDOM, VARCLICK_NONE, NULL);
}

void newdata_callback (gpointer data, guint pd_code, GtkWidget *widget) 
{
    int resp, n = 50;

    if (dataset_locked()) {
	return;
    }

    resp = spin_dialog (_("gretl: create data set"), NULL, &n, 
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

    data_structure_wizard(NULL, 1, NULL);
}

void xcorrgm_callback (gpointer p, guint v, GtkWidget *w)
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

void do_nistcheck (gpointer p, guint v, GtkWidget *w)
{
    void *handle;
    int (*run_nist_tests)(const char *, const char *, int);
    gchar *fname;
    
    run_nist_tests = gui_get_plugin_function("run_nist_tests", 
					     &handle);
    if (run_nist_tests == NULL) {
	return;
    }

    fname = g_strdup_printf("%snist.out", paths.userdir);

    (*run_nist_tests)(paths.datadir, fname, (int) v);

    close_plugin(handle);

    view_file(fname, 0, 1, 78, 400, VIEW_CODEBOOK);
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
    
    email_file(fullname, paths.userdir);
    close_plugin(handle);
}

#endif
