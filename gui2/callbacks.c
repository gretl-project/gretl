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

#ifdef OLD_GTK
# include <gtkextra/gtkiconfilesel.h>
#else
# include "treeutils.h"
#endif

extern void do_samplebool (GtkWidget *widget, dialog_t *ddata);

/* these live in dialogs.c */
extern GtkWidget *active_edit_id;
extern GtkWidget *active_edit_name;
extern GtkWidget *active_edit_text;

/* ......................................................... */

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

/* ........................................................... */

void open_data (gpointer data, guint code, GtkWidget *widget)
{
    switch (code) {
    case OPEN_DATA:
    case APPEND_DATA:
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
    } else if (action == OPEN_SESSION) {
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
    case EXPORT_DAT:
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

void print_report (gpointer data, guint unused, GtkWidget *widget)
{
    PRN *prn;

    if (bufopen(&prn)) return;

    data_report (datainfo, &paths, prn);
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

void var_resid_callback (gpointer data, guint eqnum, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data; 
    GRETL_VAR *var = (GRETL_VAR *) mydata->data;

    add_var_resid(var, eqnum);
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

void model_test_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[36], query[MAXLABEL], defstr[MAXLEN];
    char startdate[OBSLEN], enddate[OBSLEN];
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

    edit_dialog(title, query, defstr, okfunc, mydata, 
		action, varclick);   
}

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
    case SETOBS:
	strcpy(title, _("gretl: set data frequency"));
	strcpy(query, _("Enter integer frequency and\n"
	       "starting observation string:"));
	okfunc = do_setobs;
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

#ifndef OLD_GTK
    if (g_object_get_data(G_OBJECT(vwin->dialog), "text_out")) {
	u = SAVE_OUTPUT;
    } 
#else
    if (gtk_object_get_data(GTK_OBJECT(vwin->dialog), "text_out")) {
	u = SAVE_OUTPUT;
    } 
#endif
    else {
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

/* ........................................................... */

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

/* ........................................................... */

static void fix_obsstr (char *str)
{
    char pt = get_local_decpoint();
    char *p;

    p = strchr(str, ':');
    if (p != NULL) {
	*p = pt;
    }

    if (pt == ',' && (p = strchr(str, '.'))) {
	*p = pt;
    }

    if (pt == '.' && (p = strchr(str, ','))) {
	*p = pt;
    }
}

static void prep_spreadsheet (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *buf;
    char dataspec[32];
    char *test, stobs[OBSLEN], endobs[OBSLEN], firstvar[VNAMELEN];
    double sd0, ed0;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    strncpy(dataspec, buf, 31);

    /* check validity of dataspec */
    if (sscanf(dataspec, "%8s %8s %8s", stobs, endobs, firstvar) != 3) {
	errbox(_("Insufficient dataset information supplied"));
	return;
    }

    /* daily data: special */
    if (datainfo->pd == 5 || datainfo->pd == 7) {
	int err = 0;
	sd0 = (double) get_epoch_day(stobs); 
	ed0 = (double) get_epoch_day(endobs);

	if (sd0 < 0) {
	    err = 1;
	    sprintf(errtext, _("Invalid starting observation '%s'"), stobs);
	}
	if (!err && ed0 < 0) {
	    err = 1;
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	}
	if (err) {
	    errbox(errtext);
	    return;
	}
    } else { /* not daily data */
	fix_obsstr(stobs);
	fix_obsstr(endobs);

	sd0 = strtod(stobs, &test);
	if (!strcmp(stobs, test) || *test != 0 || sd0 < 0) {
	    sprintf(errtext, _("Invalid starting observation '%s'"), stobs);
	    errbox(errtext);
	    return;
	}

	ed0 = strtod(endobs, &test);
	if (!strcmp(endobs, test) || *test != 0 || ed0 < 0) {
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	    errbox(errtext);
	    return;
	}
    }

    if (sd0 > ed0) {
	sprintf(errtext, _("Empty data range '%s - %s'"), stobs, endobs);
	errbox(errtext);
	return;
    }

    colonize_obs(stobs);
    colonize_obs(endobs);

    if (datainfo->pd == 999) { /* panel */
	char unit[8], period[8];

	/* try to infer structure from ending obs */
	if (sscanf(endobs, "%[^:]:%s", unit, period) == 2) { 
	    datainfo->pd = atoi(period);
	    fprintf(stderr, I_("Setting data frequency = %d\n"), datainfo->pd);
	} else {
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	    errbox(errtext);
	    return;	    
	}
    }    

    if (datainfo->pd == 1) {
	size_t i, n;
	
	n = strlen(stobs);
	for (i=0; i<n; i++) {
	    if (!isdigit((unsigned char) stobs[i])) {
		sprintf(errtext, _("Invalid starting observation '%s'\n"
				   "for data frequency 1"), stobs);
		errbox(errtext);
		return;
	    }
	}
	n = strlen(endobs);
	for (i=0; i<n; i++) {
	    if (!isdigit((unsigned char) endobs[i])) {
		sprintf(errtext, _("Invalid ending observation '%s'\n"
				   "for data frequency 1"), endobs);
		errbox(errtext);
		return;
	    }
	}	
    } 
    else if (datainfo->pd != 5 && datainfo->pd != 7) { 
	char year[8], subper[8];

	if (sscanf(stobs, "%[^:]:%s", year, subper) != 2 ||
	    strlen(year) > 4 || atoi(subper) > datainfo->pd ||
	    (datainfo->pd < 10 && strlen(subper) != 1) ||
	    (datainfo->pd >= 10 && strlen(subper) != 2)) {
	    sprintf(errtext, _("Invalid starting observation '%s'\n"
			       "for data frequency %d"), stobs, datainfo->pd);
	    errbox(errtext);
	    return;
	}
	if (sscanf(endobs, "%[^:]:%s", year, subper) != 2 ||
	    strlen(year) > 4 || atoi(subper) > datainfo->pd ||
	    (datainfo->pd < 10 && strlen(subper) != 1) ||
	    (datainfo->pd >= 10 && strlen(subper) != 2)) {
	    sprintf(errtext, _("Invalid ending observation '%s'\n"
			       "for data frequency %d"), endobs, datainfo->pd);
	    errbox(errtext);
	    return;
	}	    
    }

    close_dialog(ddata);

    strcpy(datainfo->stobs, stobs);
    strcpy(datainfo->endobs, endobs);
    datainfo->sd0 = sd0;
    datainfo->n = -1;
    datainfo->n = dateton(datainfo->endobs, datainfo) + 1; 

    if (datainfo->n <= 0) {
	errbox("Got zero-length data series");
	return;
    }

    datainfo->v = 2;
    start_new_Z(&Z, datainfo, 0);
    datainfo->markers = 0;

    strcpy(datainfo->varname[1], firstvar);

    show_spreadsheet(datainfo);
}

/* ........................................................... */

void newdata_callback (gpointer data, guint pd_code, GtkWidget *widget) 
{
    windata_t *wdata = NULL;
    gchar *obsstr = NULL;

    if (pd_code == 0) {
	datainfo->time_series = 0;
	datainfo->pd = 1;
    } else {
	datainfo->time_series = TIME_SERIES;
	datainfo->pd = pd_code;
    }

    switch (pd_code) {
    case 0:
	datainfo->pd = 1;
	obsstr = g_strdup_printf("1 50 %s", _("newvar"));
	break;
    case 1:
	obsstr = g_strdup_printf("1950 2001 %s", _("newvar"));
	break;
    case 4:
	obsstr = g_strdup_printf("1950:1 2001:4 %s", _("newvar"));
	break;
    case 5:
	obsstr = g_strdup_printf("99/01/18 01/03/31 %s", _("newvar"));
	break;
    case 7:
	obsstr = g_strdup_printf("99/01/18 01/03/31 %s", _("newvar"));
	break;
    case 12:
	obsstr = g_strdup_printf("1950:01 2001:12 %s", _("newvar"));
	break;
    case 24:
	obsstr = g_strdup_printf("0:01 0:24 %s", _("newvar"));
	break;
    case 52:
	obsstr = g_strdup_printf("1990:01 2001:52 %s", _("newvar"));
	break;
    }

    edit_dialog (_("gretl: create data set"), 
		 _("Enter start and end obs for new data set\n"
		   "and name of first var to add:"), 
		 obsstr, 
		 prep_spreadsheet, wdata, 
		 0, 0);

    g_free(obsstr);
}

#if 0
void start_panel_callback (gpointer data, guint u, GtkWidget *widget) 
{
    windata_t *wdata = NULL;

    datainfo->pd = 999;

    edit_dialog (_("gretl: create panel data set"), 
		 _("Enter starting and ending observations and\n"
		   "the name of the first variable to add.\n"
		   "The example below is suitable for 20 units\n"
		   "observed over 10 periods"), 
		 "1:01 10:20 newvar", 
		 prep_spreadsheet, wdata, 
		 0, 0);
}
#endif

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

