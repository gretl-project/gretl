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
#include <gtkextra/gtkiconfilesel.h>

extern void browser_open_data (GtkWidget *w, gpointer data);
extern void browser_open_ps (GtkWidget *w, gpointer data);
extern void open_db_clist (GtkWidget *w, gpointer data);
extern void open_remote_clist (GtkWidget *w, gpointer data);
extern void do_samplebool (GtkWidget *widget, dialog_t *ddata);

GtkWidget *active_edit_id = NULL;
GtkWidget *active_edit_name = NULL;

/* ......................................................... */

void selectrow (GtkCList *clist, gint row, gint column, 
	        GdkEventButton *event, gpointer data) 
{
    gchar *numstr, *edttext, addvar[9];
    windata_t *mydata = (windata_t *) data;

    if (mydata == mdata) { /* main window */
	gtk_clist_get_text(GTK_CLIST(clist), row, 0, &numstr);
	mydata->active_var = atoi(numstr);
    } else {
	mydata->active_var = row;
    }

    if (active_edit_id) {
	edttext = gtk_entry_get_text (GTK_ENTRY (active_edit_id));
	if (strlen(edttext)) sprintf(addvar, " %d", mydata->active_var);
	else sprintf(addvar, "%d", mydata->active_var);
	gtk_entry_append_text(GTK_ENTRY (active_edit_id), addvar);
    }
    else if (active_edit_name) {
	edttext = gtk_entry_get_text (GTK_ENTRY (active_edit_name));
	gtk_entry_append_text(GTK_ENTRY (active_edit_name), 
			      datainfo->varname[mydata->active_var]);
	gtk_entry_append_text(GTK_ENTRY (active_edit_name), " ");
    }

    /* response to double-click */
    if (event != NULL && event->type == GDK_2BUTTON_PRESS 
	&& event->button == 1) {
	switch (mydata->role) {
	case MAINWIN:
	    display_var();
	    break;
	case RAMU_DATA:
	case GREENE_DATA:
	case JW_DATA:
	case PWT_DATA:
	    browser_open_data(NULL, mydata);
	    break;
	case RAMU_PS:
	case GREENE_PS:
	case PWT_PS:
	    browser_open_ps(NULL, mydata);
	    break;
	case NATIVE_DB:
	case RATS_DB:	    
	    open_db_clist(NULL, mydata); 
	    break;
	case REMOTE_DB:
	    open_remote_clist(NULL, mydata);
	    break;
	case NATIVE_SERIES:
	case RATS_SERIES:
	case REMOTE_SERIES:
	    gui_get_series(mydata, 0, NULL);
	    break;
	default:
	    break;
	}
    }
}

/* ........................................................... */

void open_data (gpointer data, guint code, GtkWidget *widget)
{
    switch (code) {
    case OPEN_DATA:
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
    if (action == OPEN_SCRIPT)
	file_selector(_("Open script file"), action, NULL);
    else if (action == OPEN_SESSION)
	file_selector(_("Open session file"), action, NULL);
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
    if (data_status & BOOK_DATA)
	errbox(_("You don't have permission to do this"));
    else 
	edit_buffer(&datainfo->descrip, 80, 400, _("gretl: edit data info"),
		    EDIT_HEADER);
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
    selection_dialog (_("gretl: specify model"), _("Estimate"), 
		      do_model, model_code);
}

/* ........................................................... */

#ifdef ENABLE_GMP
void mp_ols_callback (gpointer data, guint model_code, GtkWidget *widget)
{
    selection_dialog (_("gretl: specify model"), _("Estimate"), 
		      do_mp_ols, model_code);
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
    guint varclick = 0;
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
	varclick = 2;
	break;
    }

    edit_dialog(title, query, defstr, 1,
		_(" Apply "), okfunc, mydata, 
		_(" Cancel "), NULL, NULL, action, varclick);   
}

/* ........................................................... */

void add_omit_callback (gpointer data, guint action, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;

    simple_selection ("gretl: model tests", "Apply", do_add_omit, action, vwin);
}

/* ........................................................... */

void selector_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[36];
    windata_t *vwin = (windata_t *) data;
    void (*okfunc)() = NULL;

    if (action == COINT) {
	selection_dialog (_("gretl: cointegration test"), _("Estimate"), 
			  do_coint, action);
	return;
    }

    if (action == GR_XY || action == GR_IMP || action == GR_DUMMY
	|| action == SCATTERS) {
	switch (action) {
	case GR_XY:
	case GR_IMP:
	    okfunc = do_graph_from_selector;
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
	selection_dialog (_("gretl: define graph"), _("Graph"), 
			  okfunc, action);
	return;
    }

    switch (action) {
    case ADD:
    case OMIT:
	strcpy(title, _("gretl: model tests"));
	okfunc = do_add_omit;
	break;
    case PRINT:
	strcpy(title, _("gretl: display vars"));
	okfunc = display_selected;
	break;
    case GR_BOX: case GR_NBOX:
	strcpy(title, _("gretl: boxplots"));
	okfunc = do_box_graph;
	break;
    case GR_PLOT:
	strcpy(title, _("gretl: time-series plot"));
	okfunc = do_graph_from_selector;
	break;
    }

    simple_selection (title, _("Apply"), okfunc, action, vwin);
}

/* ........................................................... */

void gretl_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[64], query[MAXLABEL], defstr[MAXLEN];
    char startdate[9], enddate[9];
    void (*okfunc)() = NULL;
    int v;
    guint varclick = 0;
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
	strcpy(defstr, datainfo->varname[mdata->active_var]);
	okfunc = do_sampledum;
	varclick = 2;
	break;
    case SMPLBOOL:
	strcpy(title, _("gretl: restrict sample"));
	strcpy(query, _("Enter boolean condition for selecting cases:"));
	okfunc = do_samplebool;
	varclick = 2;
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
	varclick = 1;
	break;
    case MEANTEST:
    case MEANTEST2:
	strcpy(title, _("gretl: means test"));
	strcpy(query, _("Enter two variables by name or number:"));
	okfunc = do_dialog_cmd;
	varclick = 1;
	break;
    case VARTEST:
	strcpy(title, _("gretl: variances test"));
	strcpy(query, _("Enter two variables by name or number:"));
	okfunc = do_dialog_cmd;
	varclick = 1;
	break;
    case GENR:
	strcpy(title, _("gretl: add var"));
	strcpy(query, _("Enter formula for new variable:"));
	okfunc = do_genr;
	varclick = 2;
	break;
    case RENAME:
	strcpy(title, _("gretl: rename var"));
	strcpy(query, _("Enter new name for variable\n(max. 8 characters):"));
	strcpy(defstr, datainfo->varname[mdata->active_var]);
	okfunc = do_rename_var;
	break;
    case RELABEL:
	v = mdata->active_var;
	strcpy(title, _("gretl: edit label"));
	sprintf(query, _("Edit label for variable number %d (%s):"),
		v, datainfo->varname[v]);
	if (strlen(datainfo->label[v]) > 2) 
	    strcpy(defstr, datainfo->label[v]);
	okfunc = do_edit_label;
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
	varclick = 2;
	break;	
    default:
	errbox("Bug: unrecognized action code in gretl_callback");
	return;
    }

    edit_dialog(title, query, defstr, 1,
		_(" Apply "), okfunc, mydata, 
		_(" Cancel "), NULL, NULL, action, varclick);   
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

void file_save_callback (GtkWidget *w, gpointer data)
{
    guint u = 0;
    windata_t *mydata = (windata_t *) data;

    switch (mydata->role) {
    case EDIT_SCRIPT:
    case VIEW_SCRIPT:
	u = SAVE_SCRIPT;
	break;
    case SCRIPT_OUT:
	u = SAVE_OUTPUT;
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

    file_save(data, u, w);
}
