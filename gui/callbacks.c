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
#ifndef G_OS_WIN32
# include <gtkextra/gtkiconfilesel.h>
#endif

extern void browser_open_data (GtkWidget *w, gpointer data);
extern void browser_open_ps (GtkWidget *w, gpointer data);
extern void open_db_clist (GtkWidget *w, gpointer data);
extern void open_remote_clist (GtkWidget *w, gpointer data);
extern void do_samplebool (GtkWidget *widget, dialog_t *ddata);

GtkWidget *active_edit_id = NULL;
GtkWidget *active_edit_name = NULL;

char remember_dir[MAXLEN] = "";

/* ......................................................... */

void selectrow (GtkCList *clist, gint row, gint column, 
	        GdkEventButton *event, gpointer data) 
{
    gchar *numstr, *edttext, addvar[9];
    windata_t *mydata = (windata_t *) data;

    if (mydata == mdata) { /* main window */
	gtk_clist_get_text(GTK_CLIST(clist), row, 0, &numstr);
	mydata->active_var = atoi(numstr);
    } else 
	mydata->active_var = row;
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
    }
    /* response to double-click */
    if (event != NULL && event->type == GDK_2BUTTON_PRESS 
	&& event->button == 1) {
	switch (mydata->action) {
	case MAINWIN:
	    display_var();
	    break;
	case RAMU_DATA:
	case GREENE_DATA:
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
    char *startdir = (remember_dir[0])? remember_dir : paths.userdir;

    switch (code) {
    case OPEN_DATA:
	file_selector("Open data file", startdir, code, NULL);
	break;
    case OPEN_CSV:
	file_selector("Open CSV file", startdir, code, NULL);
	break;
    case OPEN_BOX:
	file_selector("Open BOX file", startdir, code, NULL);
	break;
    }
}

/* ........................................................... */

void open_script (gpointer data, guint action, GtkWidget *widget)
{
    switch (action) {
    case 0:
	file_selector("Open script file", paths.userdir, 
			   OPEN_SCRIPT, NULL);
	break;
    case 1:
	file_selector("Open session file", paths.userdir, 
			   OPEN_SESSION, NULL);
	break;
    case 2:
	file_selector("Open session file", paths.scriptdir, 
			   OPEN_SESSION, NULL);
	break;    
    }
}

/* ........................................................... */

void file_save (gpointer data, guint file_code, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data;
    char *startdir = (remember_dir[0])? remember_dir : paths.userdir;    

    switch (file_code) {
    case SAVE_OUTPUT:
	file_selector("Save output file", startdir, 
		      SAVE_OUTPUT, mydata->w);
	break;
    case SAVE_CONSOLE:
	file_selector("Save console output", startdir, 
		      SAVE_CONSOLE, mydata->w);
	break;
    case SAVE_CMDS: 
	file_selector("Save command log", startdir, 
		      SAVE_CMDS, mydata->w);
	break;
    case SAVE_SCRIPT:
	file_selector("Save command script", startdir, 
		      SAVE_SCRIPT, mydata->w);
	break;
    case SAVE_DATA:
	if (!storevars_dialog(STORE)) {
	    file_selector("Save data file", startdir,
			  SAVE_DATA, NULL);
	}
	break;
    case SAVE_GZDATA:
    case SAVE_BIN1:
    case SAVE_BIN2:
	if (!storevars_dialog(STORE))
	    file_selector("Save data file", startdir, 
			  file_code, NULL);
	break;
    case EXPORT_CSV:
	if (!storevars_dialog(EXPORT)) 
	    file_selector("Save CSV data file", startdir, 
			  file_code, NULL);
	break;
    case EXPORT_R:
    case EXPORT_R_ALT:
	if (!storevars_dialog(EXPORT)) 
	    file_selector("Save R data file", startdir, 
			  file_code, NULL);
	break;
    case EXPORT_OCTAVE:
	if (!storevars_dialog(EXPORT)) 
	    file_selector("Save octave data file", startdir, 
			  file_code, NULL);
	break;
    case SAVE_TEX_TAB:
    case SAVE_TEX_EQ:
	file_selector("Save LaTeX file", startdir, 
		      file_code, mydata->data);
	break;
    case SAVE_HTML:
	file_selector("Save HTML file", startdir, 
		      file_code, mydata->data);
	break;
    case SAVE_MODEL:
	file_selector("Save model output", startdir, 
		      file_code, mydata->w);
	break;
    case SAVE_GP_CMDS:
	file_selector("Save gnuplot commands", startdir, 
		      file_code, mydata->w);
	break;
    default:
	dummy_call();
    }
}

/* ........................................................... */

void dummy_call (void)
{
    errbox("Sorry, this item not yet implemented!");
}

/* ........................................................... */

void open_info (gpointer data, guint edit, GtkWidget *widget)
{
    PRN *prn;
    gint err;

    if (bufopen(&prn)) return;

    err = get_info(paths.hdrfile, prn);

    if (err) {
	gretl_print_destroy(prn);
	if (err == 1) {
	    errbox("Couldn't read data header file");
	} else if (err == 2) { /* default blank header */
	    FILE *fp;

	    fp = fopen(paths.hdrfile, "a");
	    if (fp == NULL) {
		infobox("No info in data header file");
	    } else {
		fclose(fp);
		if (!yes_no_dialog("gretl: add info", 
				  "The data header file contains no\n"
				  "informative comments.  Would you like\n"
				  "to add some now?", 0)) {
				     edit_header (NULL, 1, NULL);
		}
	    }
	}
	return;
    }
    view_buffer(prn, 80, 300, "gretl: data info", INFO, NULL);
}

/* ........................................................... */

void edit_header (gpointer data, guint save, GtkWidget *widget)
{
    view_file(paths.hdrfile, 1, 0, 80, 300, "gretl: edit header", 
	      NULL);
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
    char tempstr[MAXLEN], modelstr[MAXLEN], listnum[6];
    int i, l0 = 0, *oldlist = NULL;

    if (default_list != NULL) {
	copylist(&oldlist, default_list);
	l0 = oldlist[0];
    }

    switch (model_code) {
    case OLS:
    case HCCM:
    case HSK:
    case CORC:
    case HILU:
    case LOGIT:
    case PROBIT:
    case POOLED:
	sprintf(tempstr, "Enter specification for %s model:\n"
		"(depvar indepvars)", commands[model_code]);
	break;
    case WLS:
	strcpy(tempstr, "Enter specification for WLS model:\n"
	       "(wtvar depvar indepvars)");
	break;
    case TSLS:
	strcpy(tempstr, "Enter specification for TSLS model:\n"
	       "<varlist1> ; <varlist2>\n"
	       "You will probably want to consult the "
	       "help on tsls first");
	break;
    case AR:
	strcpy(tempstr, "Enter specification for AR model:\n"
	       "<laglist> ; <varlist>\n"
	       "You will probably want to consult the "
	       "help on ar first");
	break;
    case VAR:
	strcpy(tempstr, "Enter specification for VAR:\n"
	       "(lag_order depvar indepvars)");
	break;
    }

    if (oldlist != NULL) {
	sprintf(modelstr, "%d ", oldlist[1]);
	for (i=2; i<=l0; i++) {
	    sprintf(listnum, "%d ", oldlist[i]);
	    strcat(modelstr, listnum);
	}
    } else modelstr[0] = '\0';

    edit_dialog ("gretl: define model", tempstr,
		 modelstr, 1,
		 "Estimate", do_model, NULL, 
		 " Cancel ", NULL, NULL, model_code, 1);
    if (oldlist != NULL) free(oldlist);
}

/* ........................................................... */

void gretl_callback (gpointer data, guint action, GtkWidget *widget)
{
    char title[36], query[MAXLABEL], defstr[MAXLEN];
    char startdate[9], enddate[9];
    void (*okfunc)() = NULL;
    int v;
    guint varclick = 0;
    windata_t *mydata = (windata_t *) data;

    defstr[0] = '\0';

    switch (action) {
    case SMPL:
	strcpy(title, "gretl: set sample");
	sprintf(query, "New starting obs (min = %s)\n"
		"and ending obs (max = %s):", 
		datainfo->stobs, datainfo->endobs);
	sprintf(defstr, "%s %s", datainfo->stobs, datainfo->endobs);
	okfunc = change_sample;
	break;
    case SMPLDUM:
	strcpy(title, "gretl: define sample");
	strcpy(query, "Name of dummy variable to use:");
	strcpy(defstr, datainfo->varname[mdata->active_var]);
	okfunc = do_sampledum;
	varclick = 2;
	break;
    case SMPLBOOL:
	strcpy(title, "gretl: restrict sample");
	strcpy(query, "Enter boolean condition for selecting cases:");
	okfunc = do_samplebool;
	varclick = 2;
	break;
    case SETOBS:
	strcpy(title, "gretl: set data frequency");
	strcpy(query, "Enter integer frequency and\n"
	       "starting observation string:");
	okfunc = do_setobs;
	break;
    case SEED:
	strcpy(title, "gretl: random variables");
	strcpy(query, "Enter integer seed for\n"
	       "pseudo-random number generator:");
	okfunc = do_seed;
	break; 
    case SIM:
	if (mdata->active_var < 0) return;
	if (mdata->active_var == 0) {
	    errbox("You can't simulate with the constant");
	    return;
	}
	ntodate(startdate, datainfo->t1 + 1, datainfo);
	ntodate(enddate, datainfo->n - 1, datainfo);
	strcpy(title, "gretl: simulate variable");
	strcpy(query, "Enter spec. for simulation:\n"
	       "(start end parameters)\n"
	       "You will probably want to consult the "
	       "help on sim first");
	sprintf(defstr, "%s %s %s", startdate, enddate, 
		datainfo->varname[mdata->active_var]);
	okfunc = do_sim;
	break; 
    case NULLDATA:
	strcpy(title, "gretl: simulation data");
	strcpy(query, "Series length for simulation data set:");
	strcpy(defstr, "100");
	okfunc = do_simdata;
	break;         
    case ADF:
	strcpy(title, "gretl: ADF test");
	strcpy(query, "Lag order for ADF test:");
	strcpy(defstr, "1");
	okfunc = do_dialog_cmd;
	break;
    case ARCH:
	strcpy(title, "gretl: ARCH test");
	strcpy(query, "Lag order for ARCH test:");
	strcpy(defstr, "1");
	okfunc = do_arch;
	break;
    case COINT:
	strcpy(title, "gretl: cointegration test");
	strcpy(query, "Enter spec. for cointegration test:\n"
	       "(lag_order depvar indepvars)");
	okfunc = do_dialog_cmd;
	varclick = 1;
	break;
    case SPEARMAN:
	strcpy(title, "gretl: rank correlation");
	strcpy(query, "Enter two variables by name or number:");
	okfunc = do_dialog_cmd;
	varclick = 1;
	break;
    case MEANTEST:
    case MEANTEST2:
	strcpy(title, "gretl: means test");
	if (action == MEANTEST2) 
	    strcat(title, " (unequal var)");
	strcpy(query, "Enter two variables by name or number:");
	okfunc = do_dialog_cmd;
	varclick = 1;
	break;
    case VARTEST:
	strcpy(title, "gretl: variances test");
	strcpy(query, "Enter two variables by name or number:");
	okfunc = do_dialog_cmd;
	varclick = 1;
	break;
    case CHOW:
	ntodate(startdate, datainfo->t1, datainfo);
	ntodate(enddate, datainfo->t2, datainfo);
	strcpy(title, "gretl: Chow test");
	sprintf(query, "Enter observation at which\n"
		"to split the sample\n"
		"(between %s and %s):", startdate, enddate);
	okfunc = do_chow;
	break;
    case OMIT:
	strcpy(title, "gretl: omit vars");
	strcpy(query, "Names (or numbers) of variables to omit:");
	okfunc = do_add_omit;
	varclick = 1;
	break;	
    case ADD:
	strcpy(title, "gretl: add vars");
	strcpy(query, "Names (or numbers) of variables to add:");
	okfunc = do_add_omit;
	varclick = 1;
	break;
    case FCAST: 
	strcpy(title, "gretl: forecast");
	sprintf(query, "Starting obs (min = %s)\n"
		"and ending obs (max = %s)?", 
		datainfo->stobs, datainfo->endobs);
	sprintf(defstr, "%s %s", datainfo->stobs, datainfo->endobs);
	okfunc = do_forecast;
	break;
    case PRINT:
	strcpy(title, "gretl: display vars");
	strcpy(query, "Enter variable names or numbers:");
	okfunc = display_selected;
	varclick = 1;
	break;
    case GENR:
	strcpy(title, "gretl: add var");
	strcpy(query, "Enter formula for new variable:");
	okfunc = do_genr;
	varclick = 2;
	break;
    case MODEL_GENR:
	strcpy(title, "gretl: add var");
	strcpy(query, "Enter formula for new variable:");
	okfunc = do_model_genr;
	varclick = 2;
	break;
    case RENAME:
	strcpy(title, "gretl: rename var");
	strcpy(query, "Enter new name for variable\n(max. 8 characters):");
	strcpy(defstr, datainfo->varname[mdata->active_var]);
	okfunc = do_rename_var;
	break;
    case RELABEL:
	v = mdata->active_var;
	strcpy(title, "gretl: edit label");
	sprintf(query, "Edit label for variable number %d (%s):",
		v, datainfo->varname[v]);
	if (strlen(datainfo->label[v]) > 2) 
	    strcpy(defstr, datainfo->label[v]);
	okfunc = do_edit_label;
	break;
    case MARKERS:
	strcpy(title, "gretl: add markers");
	strcpy(query, "Supply full path to file with markers:");
	strcpy(defstr, paths.userdir);
	okfunc = do_add_markers; 
	break;
    case CORRGM:
	strcpy(title, "gretl: correlogram");
	strcpy(query, "Max lag length?\n(0 for automatic):");
	strcpy(defstr, "0");
	okfunc = do_dialog_cmd;
	break;
    }

    edit_dialog(title, query, defstr, 1,
		" Apply ", okfunc, mydata, 
		" Cancel ", NULL, NULL, action, varclick);   
}

