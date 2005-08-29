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
#include "ssheet.h"
#include "textbuf.h"
#include "textutil.h"
#include "boxplots.h"
#include "dlgutils.h"

#ifdef OLD_GTK
# include <gtkextra/gtkiconfilesel.h>
#else
# include "treeutils.h"
#endif

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
    case NATIVE_DB:
    case RATS_DB:	    
	open_db_list(NULL, win); 
	break;
    case REMOTE_DB:
	open_remote_db_list(NULL, win);
	break;
    case NATIVE_SERIES:
    case RATS_SERIES:
    case REMOTE_SERIES:
	gui_get_series(win, 0, NULL);
	break;
    default:
	break;
    }
}

#ifndef OLD_GTK

void listbox_select_row (GtkTreeSelection *selection, gpointer data)
{
    windata_t *win = (windata_t *) data;
    GtkTreeIter iter;
    GtkTreeModel *model;
    GtkTreePath *path;
    gint row;

    if (!gtk_tree_selection_get_selected (selection, &model, &iter))
	return;

    path = gtk_tree_model_get_path (model, &iter);
    row = tree_path_get_row_number (path);
    win->active_var = row;
    gtk_tree_path_free (path);
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

#else /* now comes an old gtk function */

void selectrow (GtkCList *clist, gint row, gint column, 
	        GdkEventButton *event, gpointer data) 
{
    gchar *numstr, *edttext;
    windata_t *win = (windata_t *) data;

    if (win == mdata) { /* main window */
	gtk_clist_get_text(clist, row, 0, &numstr);
	win->active_var = atoi(numstr);
    } else {
	win->active_var = row;
    }

    if (active_edit_id != NULL) {
	gchar addvar[VNAMELEN];

	edttext = gtk_entry_get_text(GTK_ENTRY(active_edit_id));
	if (*edttext != '\0') {
	    sprintf(addvar, " %d", win->active_var);
	} else {
	    sprintf(addvar, "%d", win->active_var);
	}
	gtk_entry_append_text(GTK_ENTRY(active_edit_id), addvar);
    } else if (active_edit_name != NULL) {
	edttext = gtk_entry_get_text (GTK_ENTRY(active_edit_name));
	gtk_entry_append_text(GTK_ENTRY(active_edit_name), 
			      datainfo->varname[win->active_var]);
	gtk_entry_append_text(GTK_ENTRY(active_edit_name), " ");
    }

    /* response to double-click */
    if (event != NULL && event->type == GDK_2BUTTON_PRESS 
	&& event->button == 1) {
	doubleclick_action(win);
    }
}

void unselectrow (GtkCList *clist, gint row, gint column, 
		  GdkEventButton *event, gpointer data) 
{
    windata_t *win = (windata_t *) data;

    if (win != mdata) { /* main window */
	return;
    } else {
	gchar *numstr;

	gtk_clist_get_text(clist, row, 0, &numstr);
	if (win->active_var == atoi(numstr)) {
	    win->active_var = clist->focus_row;
	}
    }
}

#endif /* old versus new gtk */

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
	file_selector(_("Open data file"), code, NULL);
	break;
    case OPEN_CSV:
    case APPEND_CSV:
	delimiter_dialog(NULL);
	file_selector(_("Open CSV file"), code, NULL);
	break;
    case OPEN_OCTAVE:
    case APPEND_OCTAVE:
	file_selector(_("Open Octave file"), code, NULL);
	break;
    case OPEN_BOX:
	file_selector(_("Open BOX file"), code, NULL);
	break;
    case OPEN_GNUMERIC:
    case APPEND_GNUMERIC:
	file_selector(_("Open Gnumeric file"), code, NULL);
	break;
    case OPEN_EXCEL:
    case APPEND_EXCEL:
	file_selector(_("Open Excel file"), code, NULL);
	break;
    case OPEN_WF1:
    case APPEND_WF1:
	file_selector(_("Open Eviews workfile"), code, NULL);
	break;
    case OPEN_DTA:
    case APPEND_DTA:
	file_selector(_("Open Stata file"), code, NULL);
	break;
    case OPEN_MARKERS:
	file_selector(_("gretl: add markers"), code, NULL);
	break;
    default:
	errbox("Unrecognized data code");
	break;
    }
}

void open_script (gpointer data, guint action, GtkWidget *widget)
{
    if (action == OPEN_SCRIPT) {
	file_selector(_("Open script file"), action, NULL);
    } else if (action == OPEN_SESSION) {
	file_selector(_("Open session file"), action, NULL);
    }
}

void file_save (gpointer data, guint file_code, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    gretlopt opt = OPT_NONE;
    gpointer p = NULL;

    switch (file_code) {
    case SAVE_OUTPUT:
	file_selector(_("Save output file"), SAVE_OUTPUT, vwin);
	break;
    case SAVE_CONSOLE:
	file_selector(_("Save console output"), SAVE_CONSOLE, vwin);
	break;
    case SAVE_CMDS: 
	file_selector(_("Save command log"), SAVE_CMDS, vwin);
	break;
    case SAVE_SCRIPT:
	file_selector(_("Save command script"), SAVE_SCRIPT, vwin);
	break;
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case SAVE_GZDATA:
    case SAVE_BIN1:
    case SAVE_BIN2:
    case SAVE_DBDATA:
	data_save_selection_wrapper(file_code, NULL);
	break;
    case EXPORT_CSV:
	delimiter_dialog(&opt);
	p = GINT_TO_POINTER(opt);
    case EXPORT_R:
    case EXPORT_R_ALT:
    case EXPORT_OCTAVE:
    case EXPORT_DAT:
	data_save_selection_wrapper(file_code, p);
	break;
    case SAVE_TEX:
	file_selector(_("Save LaTeX file"), file_code, vwin->data);
	break;
    case SAVE_MODEL:
	file_selector(_("Save model output"), file_code, vwin);
	break;
    case SAVE_GP_CMDS:
	file_selector(_("Save gnuplot commands"), file_code, vwin);
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
    windata_t *mydata = (windata_t *) data; 
    GRETL_VAR *var = (GRETL_VAR *) mydata->data;

    add_var_resid(var, eqnum);
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

#ifdef ENABLE_GMP
void mp_ols_callback (gpointer data, guint model_code, GtkWidget *widget)
{
    selection_dialog(_("gretl: specify model"), do_mp_ols, model_code, 0);
}
#endif /* ENABLE_GMP */

void model_genr_callback (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data;

    edit_dialog(_("gretl: add var"), _("Enter formula for new variable:"),
		"", do_model_genr, mydata, 
		MODEL_GENR, VARCLICK_INSERT_NAME);   
}

void selector_callback (gpointer data, guint action, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    char title[64];

    if (action == ADD || action == OMIT || action == COEFFSUM) {
	set_window_busy(vwin);
    }

    strcpy(title, "gretl: ");

    if (action == COINT || action == COINT2) {
	selection_dialog(_("gretl: cointegration test"), do_coint, action, 0);
    } else if (action == VAR || action == VECM) {
	selection_dialog((action == VAR)? _("gretl: VAR") : _("gretl: VECM"),
			 do_vector_model, action, 0);
    } else if (action == GR_XY || action == GR_IMP || action == GR_DUMMY
	       || action == SCATTERS || action == GR_3D) {
	void (*okfunc)() = NULL;

	switch (action) {
	case GR_XY:
	case GR_IMP:
	    okfunc = do_graph_from_selector;
	    break;
	case GR_3D:
	    okfunc = do_splot_from_selector;
	    break;
	case GR_DUMMY:
	    okfunc = do_dummy_graph;
	    break;
	case SCATTERS:
	    okfunc = do_scatters;
	    break;
	default:
	    return;
	}
	selection_dialog(_("gretl: define graph"), okfunc, action, 0);
    } else if (action == ADD || action == OMIT) {
	simple_selection(_("gretl: model tests"), do_add_omit, action, vwin);
    } else if (action == COEFFSUM) {
	simple_selection(_("gretl: model tests"), do_coeff_sum, action, vwin);
    } else if (action == GR_PLOT) {
	simple_selection(_("gretl: model tests"), do_graph_from_selector, action, vwin);
    } else if (action == SPEARMAN) {
	strcat(title, _("rank correlation"));
	simple_selection(title, do_spearman, action, vwin);
    } else if (action == MEANTEST || action == MEANTEST2) {
	strcpy(title, _("gretl: means test"));
	simple_selection(title, do_two_var_test, action, vwin);
    } else if (action == VARTEST) {
	strcpy(title, _("gretl: variances test"));
	simple_selection(title, do_two_var_test, action, vwin);
    } else {
	errbox("selector_callback: code was not recognized");
    }
}

void gretl_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[64], query[MAXLABEL], defstr[MAXLEN];
    void (*okfunc)() = NULL;
    guint varclick = VARCLICK_NONE;
    windata_t *mydata = (windata_t *) data;

    *defstr = '\0';

    switch (action) {
    case SMPLBOOL:
	strcpy(title, _("gretl: restrict sample"));
	strcpy(query, _("Enter boolean condition for selecting cases:"));
	okfunc = do_samplebool;
	varclick = VARCLICK_INSERT_NAME;
	break;
    case SETSEED:
	strcpy(title, _("gretl: random variables"));
	strcpy(query, _("Enter integer seed for\n"
	       "pseudo-random number generator:"));
	okfunc = do_seed;
	break; 
    case NULLDATA:
	strcpy(title, _("gretl: simulation data"));
	strcpy(query, _("Series length for simulation data set:"));
	strcpy(defstr, "100");
	okfunc = do_simdata;
	break;         
    case GENR:
	strcpy(title, _("gretl: add var"));
	strcpy(query, _("Enter formula for new variable:"));
	okfunc = do_genr;
	varclick = VARCLICK_INSERT_NAME;
	break;
    case VSETMISS:
	strcpy(title, _("gretl: missing code"));
	strcpy(query, _("Enter value to be read as \"missing\":"));
	okfunc = do_variable_setmiss;
	break;
    case GSETMISS:
	strcpy(title, _("gretl: missing code"));
	strcpy(query, _("Enter value to be read as \"missing\":"));
	okfunc = do_global_setmiss;
	break;
    case GR_BOX:
    case GR_NBOX:
	strcpy(title, _("gretl: boxplots"));
	strcpy(query, _("Specify variables to plot:"));
	okfunc = do_box_graph;
	varclick = VARCLICK_INSERT_NAME;
	strcpy(defstr, get_boxplots_string());
	break;	
    case NLS:
	strcpy(title, _("gretl: nonlinear least squares"));
	strcpy(query, _("NLS: Specify function, and derivatives if possible:"));
	okfunc = do_nls_model;
	varclick = VARCLICK_INSERT_TEXT;
	break;	
    case RESTRICT:
	strcpy(title, _("gretl: linear restrictions"));
	strcpy(query, _("Specify restrictions:"));
	okfunc = do_restrict;
	break;	
    default:
	errbox("Bug: unrecognized action code in gretl_callback");
	return;
    }

    edit_dialog(title, query, defstr, okfunc, mydata, 
		action, varclick);   
}

void text_copy_callback (GtkWidget *w, gpointer data)
{
    window_copy(data, GRETL_FORMAT_SELECTION, w);
}

void text_paste_callback (GtkWidget *w, gpointer data)
{
    text_paste(data, 0, w);
}

void text_replace_callback (GtkWidget *w, gpointer data)
{
    text_replace(data, 0, w);
}

void text_undo_callback (GtkWidget *w, gpointer data)
{
    text_undo(data, 0, w);
}

void run_script_callback (GtkWidget *w, gpointer data)
{
    do_run_script(data, SCRIPT_EXEC, w);
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
	    u = SAVE_SCRIPT;
	    break;
	case VIEW_LOG:
	    u = SAVE_CMDS;
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

void add_random_callback (gpointer data, guint code, GtkWidget *widget) 
{
    if (code == GENR_UNIFORM) {
	edit_dialog (_("gretl: uniform variable"), 
		     _("Enter name for variable, and\n"
		       "minimum and maximum values:"), 
		     "unif 0 1",  
		     do_random, NULL, 
		     GENR_UNIFORM, GENR);
    } else if (code == GENR_NORMAL) {
	edit_dialog (_("gretl: normal variable"), 
		     _("Enter name, mean and standard deviation:"), 
		     "norm 0 1", 
		     do_random, NULL, 
		     GENR_NORMAL, GENR);
    }
}

static void name_first_ssheet_var (GtkWidget *widget, dialog_t *dlg) 
{
    const gchar *buf;

    buf = edit_dialog_get_text(dlg);

    if (buf == NULL || validate_varname(buf)) {
	return;
    }

    datainfo->varname[1][0] = 0;
    strncat(datainfo->varname[1], buf, VNAMELEN - 1);

    close_dialog(dlg);

    show_spreadsheet(SHEET_NEW_DATASET);
}

static int prep_spreadsheet (const char *dataspec)
{
    char stobs[OBSLEN], endobs[OBSLEN];
    int t1, t2;

    if (sscanf(dataspec, "%10s %10s", stobs, endobs) != 2) {
	errbox(_("Insufficient dataset information supplied"));
	return 1;
    }

    t1 = dateton(stobs, datainfo);
    t2 = dateton(endobs, datainfo);

    strcpy(datainfo->stobs, stobs);
    strcpy(datainfo->endobs, endobs);

    datainfo->sd0 = get_date_x(datainfo->pd, datainfo->stobs);

    datainfo->n = t2 - t1 + 1;
    datainfo->t1 = 0;
    datainfo->t2 = datainfo->n - 1;

    datainfo->v = 2;
    start_new_Z(&Z, datainfo, 0);
    datainfo->markers = 0;

    edit_dialog (_("gretl: name variable"), 
		 _("Enter name for new variable\n"
		   "(max. 8 characters)"),
		 NULL, name_first_ssheet_var, NULL, 
		 0, 0);

    return 0;
}

static void n_obs_callback (GtkWidget *w, dialog_t *dlg)
{
    const gchar *buf;
    char obsstr[32];
    int n;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    n = atoi(buf);
    if (n < 1 || n > 10000) {
	errbox(_("Invalid number of observations"));
	return;
    }

    close_dialog(dlg);

    datainfo->n = n;
    sprintf(obsstr, "1 %d", n);
    prep_spreadsheet(obsstr);
}

void newdata_callback (gpointer data, guint pd_code, GtkWidget *widget) 
{
    if (dataset_locked()) {
	return;
    }

    if (pd_code == 0) {
	/* cross-sectional dataset */
	datainfo->structure = CROSS_SECTION;
	datainfo->pd = 1;
	datainfo->sd0 = 1.0;
	strcpy(datainfo->stobs, "1");
	edit_dialog (_("gretl: create data set"), 
		     _("Number of observations:"), "50",
		     n_obs_callback, NULL, 0, 0);
    } else {
	gchar *obsstr = NULL;

	datainfo->structure = TIME_SERIES;
	datainfo->pd = pd_code;
	compute_default_ts_info(datainfo, 1);

	sample_range_dialog(&obsstr, CREATE_DATASET, NULL);

	if (obsstr != NULL) {
	    prep_spreadsheet(obsstr);
	    g_free(obsstr);
	}
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
