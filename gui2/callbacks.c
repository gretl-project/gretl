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
#include "treeutils.h"
#include "selector.h"
#include "session.h"

extern void browser_open_data (GtkWidget *w, gpointer data);
extern void browser_open_ps (GtkWidget *w, gpointer data);
extern void open_db_list (GtkWidget *w, gpointer data);
extern void open_remote_db_list (GtkWidget *w, gpointer data);
extern void do_samplebool (GtkWidget *widget, dialog_t *ddata);

/* ......................................................... */

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

/* ......................................................... */

gint listbox_double_click (GtkWidget *widget, GdkEventButton *event,
			   windata_t *win)
{
    if (event != NULL && event->type == GDK_2BUTTON_PRESS 
	&& event->button == 1) {
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
    return FALSE;
}

/* ........................................................... */

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

/* ........................................................... */

void open_data (gpointer data, guint code, GtkWidget *widget)
{
    switch (code) {
    case OPEN_DATA:
    case OPEN_ASCII:
    case APPEND_ASCII:
	file_selector(_("Open data file"), code, NULL);
	break;
    case OPEN_CSV:
    case APPEND_CSV:
	delimiter_dialog();
	file_selector(_("Open CSV file"), code, NULL);
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
    case OPEN_DES:
	 file_selector(_("Open Wooldridge file"), code, NULL);
         break;
    default:
	errbox("Unrecognized data code");
	break;
    }
}

/* ........................................................... */

void open_script (gpointer data, guint action, GtkWidget *widget)
{
    if (action == OPEN_SCRIPT) {
	file_selector(_("Open script file"), action, NULL);
    }
    else if (action == OPEN_SESSION) {
	file_selector(_("Open session file"), action, NULL);
    }
}

/* ........................................................... */

void file_save (gpointer data, guint file_code, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data;

    switch (file_code) {
    case SAVE_OUTPUT:
	file_selector(_("Save output file"), SAVE_OUTPUT, mydata);
	break;
    case SAVE_CONSOLE:
	file_selector(_("Save console output"), SAVE_CONSOLE, mydata);
	break;
    case SAVE_CMDS: 
	file_selector(_("Save command log"), SAVE_CMDS, mydata);
	break;
    case SAVE_SCRIPT:
	file_selector(_("Save command script"), SAVE_SCRIPT, mydata);
	break;
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case SAVE_GZDATA:
    case SAVE_BIN1:
    case SAVE_BIN2:
	data_save_selection_wrapper(file_code);
	break;
    case EXPORT_CSV:
	delimiter_dialog();
    case EXPORT_R:
    case EXPORT_R_ALT:
    case EXPORT_OCTAVE:
	data_save_selection_wrapper(file_code);
	break;
    case SAVE_TEX_TAB:
    case SAVE_TEX_EQ:
    case SAVE_TEX_TAB_FRAG:
    case SAVE_TEX_EQ_FRAG:
	file_selector(_("Save LaTeX file"), file_code, mydata->data);
	break;
    case SAVE_MODEL:
	file_selector(_("Save model output"), file_code, mydata);
	break;
    case SAVE_GP_CMDS:
	file_selector(_("Save gnuplot commands"), file_code, mydata);
	break;
    default:
	dummy_call();
    }
}

/* ........................................................... */

void dummy_call (void)
{
    errbox(_("Sorry, this item not yet implemented!"));
}

/* ........................................................... */

/* contortions are needed here to get around the fact that the
   output of strftime (used in print_time()) will not be UTF-8 */

/* actually, I now think the above is mistaken -- Oct 02 */

void print_report (gpointer data, guint unused, GtkWidget *widget)
{
    PRN *prn;
#ifdef BROKEN_NLS /* was ENABLE_NLS */ 
    gchar *utfbuf;
    gsize wrote;
#endif  

    if (bufopen(&prn)) return;

#ifdef BROKEN_NLS
    bind_textdomain_codeset(PACKAGE, "ISO-8859-1");
#endif
    data_report (datainfo, &paths, prn);
#ifdef BROKEN_NLS
    bind_textdomain_codeset(PACKAGE, "UTF-8");
    utfbuf = g_locale_to_utf8(prn->buf, -1, NULL, &wrote, NULL);
    if (utfbuf != NULL) {
	free(prn->buf);
	prn->buf = utfbuf;
    }
#endif    
    view_buffer(prn, 77, 400, _("gretl: data summary"), 
		DATA_REPORT, NULL);
}

/* ........................................................... */

void edit_header (gpointer data, guint unused, GtkWidget *widget)
{
    if (data_status & BOOK_DATA) {
	errbox(_("You don't have permission to do this"));
    } else { 
	edit_buffer(&datainfo->descrip, 80, 400, _("gretl: edit data info"),
		    EDIT_HEADER);
    }
}

/* ........................................................... */

void fit_resid_callback (gpointer data, guint code, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data; 
    MODEL *pmod = mydata->data;

    add_fit_resid(pmod, code, 0);
}

/* ........................................................... */

void model_stat_callback (gpointer data, guint which, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data; 
    MODEL *pmod = mydata->data;

    add_model_stat(pmod, which);
}

/* ........................................................... */

void model_callback (gpointer data, guint model_code, GtkWidget *widget) 
{
    selection_dialog (_("gretl: specify model"), do_model, model_code);
}

/* ........................................................... */

#ifdef ENABLE_GMP
void mp_ols_callback (gpointer data, guint model_code, GtkWidget *widget)
{
    selection_dialog (_("gretl: specify model"), do_mp_ols, model_code);
}
#endif /* ENABLE_GMP */

/* ........................................................... */

#ifdef notdef
static int model_dates_check (windata_t *mydata)
{
    MODEL *pmod = (MODEL *) mydata->data;

    if (pmod->smpl.t1 != datainfo->t1 ||
	pmod->smpl.t2 != datainfo->t2) {
	errbox(_("Sorry, can't do: the sample has been reset\nsince this model "
	       "was estimated"));
	return 1;
    }
    return 0;
}
#endif

/* ........................................................... */

void model_test_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[36], query[MAXLABEL], defstr[MAXLEN];
    char startdate[9], enddate[9];
    void (*okfunc)() = NULL;
    guint varclick = VARCLICK_NONE;
    windata_t *mydata = (windata_t *) data;

    *defstr = '\0';

    switch (action) {
    case ARCH:
	strcpy(title, _("gretl: ARCH test"));
	strcpy(query, _("Lag order for ARCH test:"));
	strcpy(defstr, "1");
	okfunc = do_arch;
	break;
    case LMTEST: 
	strcpy(title, _("gretl: autocorrelation"));
	strcpy(query, _("Lag order for test:"));
	if (dataset_is_panel(datainfo)) {
	    strcpy(defstr, "1");
	} else {
	    sprintf(defstr, "%d", datainfo->pd);
	}
	okfunc = do_autocorr;
	break;
    case CHOW:
	ntodate(startdate, datainfo->t1, datainfo);
	ntodate(enddate, datainfo->t2, datainfo);
	strcpy(title, _("gretl: Chow test"));
	sprintf(query, _("Enter observation at which\n"
		"to split the sample\n"
		"(between %s and %s):"), startdate, enddate);
	okfunc = do_chow;
	break;
    case FCASTERR: 
	strcpy(title, _("gretl: forecast"));
	sprintf(query, _("Starting obs (min = %s)\n"
		"and ending obs (max = %s)?"), 
		datainfo->stobs, datainfo->endobs);
	sprintf(defstr, "%s %s", datainfo->stobs, datainfo->endobs);
	okfunc = do_forecast;
	break;
    case MODEL_GENR:
	strcpy(title, _("gretl: add var"));
	strcpy(query, _("Enter formula for new variable:"));
	okfunc = do_model_genr;
	varclick = VARCLICK_INSERT_NAME;
	break;
    }

    edit_dialog(title, query, defstr, 
		okfunc, mydata, 
		action, varclick);   
}

/* ........................................................... */

#ifdef notdef
void add_omit_callback (gpointer data, guint action, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;

    simple_selection ("gretl: model tests", do_add_omit, action, vwin);
}
#endif

/* ........................................................... */

void selector_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[36];
    windata_t *vwin = (windata_t *) data;
    void (*okfunc)() = NULL;

    if (action == COINT) {
	selection_dialog (_("gretl: cointegration test"), do_coint, action);
	return;
    }

    if (action == COINT2) {
	selection_dialog (_("gretl: cointegration test"), do_coint2, action);
	return;
    }

    if (action == GR_XY || action == GR_IMP || action == GR_DUMMY
	|| action == SCATTERS || action == GR_3D) {
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
	selection_dialog (_("gretl: define graph"), okfunc, action);
	return;
    }

    if (action == ADD || action == OMIT) {
	strcpy(title, _("gretl: model tests"));
	simple_selection (title, do_add_omit, action, vwin);
	return;
    }

    if (action == COEFFSUM) {
	strcpy(title, _("gretl: model tests"));
	simple_selection (title, do_coeff_sum, action, vwin);
	return;
    }

    if (action == GR_PLOT) {
	strcpy(title, _("gretl: model tests"));
	simple_selection (title, do_graph_from_selector, action, vwin);
	return;
    }

    errbox("selector_callback: code was not recognized");
}

/* ........................................................... */

static void maybe_insert_varname (char *s)
{
    int v = mdata->active_var;

    if (datainfo->vector[v] && isdummy(Z[v], 0, datainfo->n)) {
	strcpy(s, datainfo->varname[v]);
    }
}

/* ........................................................... */

void gretl_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[64], query[MAXLABEL], defstr[MAXLEN];
    char startdate[9], enddate[9];
    void (*okfunc)() = NULL;
    guint varclick = VARCLICK_NONE;
    windata_t *mydata = (windata_t *) data;

    defstr[0] = '\0';

    switch (action) {
    case SMPL:
	strcpy(title, _("gretl: set sample"));
	sprintf(query, _("New starting obs (min = %s)\n"
		"and ending obs (max = %s):"), 
		datainfo->stobs, datainfo->endobs);
	sprintf(defstr, "%s %s", datainfo->stobs, datainfo->endobs);
	okfunc = change_sample;
	break;
    case SMPLDUM:
	strcpy(title, _("gretl: define sample"));
	strcpy(query, _("Name of dummy variable to use:"));
	maybe_insert_varname(defstr);
	okfunc = do_sampledum;
	varclick = VARCLICK_INSERT_NAME;
	break;
    case SMPLBOOL:
	strcpy(title, _("gretl: restrict sample"));
	strcpy(query, _("Enter boolean condition for selecting cases:"));
	okfunc = do_samplebool;
	varclick = VARCLICK_INSERT_NAME;
	break;
    case SETOBS:
	strcpy(title, _("gretl: set data frequency"));
	strcpy(query, _("Enter integer frequency and\n"
	       "starting observation string:"));
	okfunc = do_setobs;
	break;
    case SEED:
	strcpy(title, _("gretl: random variables"));
	strcpy(query, _("Enter integer seed for\n"
	       "pseudo-random number generator:"));
	okfunc = do_seed;
	break; 
    case SIM:
	if (mdata->active_var < 0) return;
	if (mdata->active_var == 0) {
	    errbox(_("You can't simulate with the constant"));
	    return;
	}
	ntodate(startdate, datainfo->t1 + 1, datainfo);
	ntodate(enddate, datainfo->n - 1, datainfo);
	strcpy(title, _("gretl: simulate variable"));
	strcpy(query, _("Enter spec. for simulation:\n"
	       "(start end parameters)\n"
	       "You will probably want to consult the "
	       "help on sim first"));
	sprintf(defstr, "%s %s %s", startdate, enddate, 
		datainfo->varname[mdata->active_var]);
	okfunc = do_sim;
	break; 
    case NULLDATA:
	strcpy(title, _("gretl: simulation data"));
	strcpy(query, _("Series length for simulation data set:"));
	strcpy(defstr, "100");
	okfunc = do_simdata;
	break;         
    case ADF:
	strcpy(title, _("gretl: ADF test"));
	strcpy(query, _("Lag order for ADF test:"));
	strcpy(defstr, "1");
	okfunc = do_dialog_cmd;
	break;
    case SPEARMAN:
	strcpy(title, _("gretl: rank correlation"));
	strcpy(query, _("Enter two variables by name or number:"));
	okfunc = do_dialog_cmd;
	varclick = VARCLICK_INSERT_ID;
	break;
    case MEANTEST:
    case MEANTEST2:
	strcpy(title, _("gretl: means test"));
	strcpy(query, _("Enter two variables by name or number:"));
	okfunc = do_dialog_cmd;
	varclick = VARCLICK_INSERT_ID;
	break;
    case VARTEST:
	strcpy(title, _("gretl: variances test"));
	strcpy(query, _("Enter two variables by name or number:"));
	okfunc = do_dialog_cmd;
	varclick = VARCLICK_INSERT_ID;
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
    case MARKERS:
	strcpy(title, _("gretl: add markers"));
	strcpy(query, _("Supply full path to file with markers:"));
	strcpy(defstr, paths.userdir);
	okfunc = do_add_markers; 
	break;
    case CORRGM:
	strcpy(title, _("gretl: correlogram"));
	strcpy(query, _("Max lag length?\n(0 for automatic):"));
	strcpy(defstr, "0");
	okfunc = do_dialog_cmd;
	break;
    case GR_BOX:
    case GR_NBOX:
	strcpy(title, _("gretl: boxplots"));
	strcpy(query, _("Specify variables to plot:"));
	okfunc = do_box_graph_trad;
	varclick = VARCLICK_INSERT_NAME;
	break;	
    case NLS:
	strcpy(title, _("gretl: nonlinear least squares"));
	strcpy(query, _("NLS: Specify function, and derivatives if possible:"));
	okfunc = do_nls_model;
	varclick = VARCLICK_INSERT_TEXT;
	break;	
    default:
	errbox("Bug: unrecognized action code in gretl_callback");
	return;
    }

    edit_dialog(title, query, defstr, 
		okfunc, mydata, 
		action, varclick);   
}

/* ........................................................... */

void delete_var_by_id (int id)
{
    if (dataset_drop_var(id, &Z, datainfo))
	errbox(_("Out of memory reorganizing data set"));
    else {
	refresh_data();
	if (id < datainfo->v - 1)
	    infobox(_("Take note: variables have been renumbered"));
    }
}

/* ........................................................... */

void text_copy_callback (GtkWidget *w, gpointer data)
{
    text_copy(data, COPY_SELECTION, w);
}

/* ........................................................... */

void text_paste_callback (GtkWidget *w, gpointer data)
{
    text_paste(data, 0, w);
}

/* ........................................................... */

void text_replace_callback (GtkWidget *w, gpointer data)
{
    text_replace(data, 0, w);
}

/* ........................................................... */

void text_undo_callback (GtkWidget *w, gpointer data)
{
    text_undo(data, 0, w);
}

/* ........................................................... */

void run_script_callback (GtkWidget *w, gpointer data)
{
    do_run_script(data, SCRIPT_EXEC, w);
}

/* ........................................................... */

void gp_send_callback (GtkWidget *w, gpointer data)
{
    gp_to_gnuplot(data, 0, w);
}

/* ........................................................... */

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
