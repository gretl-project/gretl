/*
 *  Copyright (c) by Allin Cottrell
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

/* fancy_dialog.c for gretl */

#include "gretl.h"
#include "fancy_dialog.h"

static void set_weight_var (gint i, new_dialog *nd)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(nd->extra_entry), vname);
    gtk_object_set_user_data(GTK_OBJECT(nd->extra_entry),
			     GINT_TO_POINTER(atoi(vnum)));
}

static void select_weight_callback (GtkWidget *w, new_dialog *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) set_weight_var, nd);
    }
}

static void set_dependent_var (gint i, new_dialog *nd)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(nd->depvar), vname);
    gtk_object_set_user_data(GTK_OBJECT(nd->depvar),
			     GINT_TO_POINTER(atoi(vnum)));
}

static void select_dependent_callback (GtkWidget *w, new_dialog *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) set_dependent_var, nd);
    }
}

static void add_independent_var (gint i, new_dialog *nd)
{
    gchar *row[2];
    gint j, rows = GTK_CLIST(nd->indepvars)->rows;
    gint already_there = 0;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &row[0]);
    for (j=0; j<rows; j++) {
	gchar *test;

	gtk_clist_get_text(GTK_CLIST(nd->indepvars), j, 0, &test);
	if (!strcmp(test, row[0])) {
	    already_there = 1; 
	    break;
	}
    }
    if (!already_there) {
	gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 1, &row[1]);
	gtk_clist_append(GTK_CLIST(nd->indepvars), row);
    }
}

static void add_independent_callback (GtkWidget *w, new_dialog *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) add_independent_var, nd);
}

static void remove_independent_var (gint i, new_dialog *nd)
{
    gtk_clist_remove(GTK_CLIST(nd->indepvars), i);
}

static void remove_independent_callback (GtkWidget *w, new_dialog *nd)
{
    GList *mylist = GTK_CLIST(nd->indepvars)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) remove_independent_var, nd);
}

static void clear_vars (GtkWidget *w, new_dialog *nd)
{
    gchar *row[2];

    gtk_clist_unselect_all(GTK_CLIST(nd->varlist));    
    gtk_entry_set_text(GTK_ENTRY(nd->depvar), ""); 
    gtk_clist_clear(GTK_CLIST(nd->indepvars));
    row[0] = "0";
    row[1] = "const";
    gtk_clist_append(GTK_CLIST(nd->indepvars), row);
}

static void construct_cmdlist (GtkWidget *w, new_dialog *nd)
{
    gint i, rows = GTK_CLIST(nd->indepvars)->rows;
    gchar numstr[6];

    nd->cmdlist[0] = 0;

    if (nd->code == WLS) {
	i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(nd->extra_entry)));
	sprintf(numstr, "%d ", i);
	strcat(nd->cmdlist, numstr);
    }
    else if (nd->code == AR) {
	gchar *lags;

	lags = gtk_entry_get_text(GTK_ENTRY(nd->extra_entry));
	strcat(nd->cmdlist, lags);
	strcat(nd->cmdlist, " ; ");
    }

    i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(nd->depvar)));
    sprintf(numstr, "%d", i);
    strcat(nd->cmdlist, numstr);
    
    if (GTK_TOGGLE_BUTTON(nd->default_check)->active) 
	*(nd->default_var) = i;

    for (i=0; i<rows; i++) {
	gchar *indep;

	gtk_clist_get_text(GTK_CLIST(nd->indepvars), i, 0, &indep);
	strcat(nd->cmdlist, " ");
	strcat(nd->cmdlist, indep);
    }
    fprintf(stderr, "cmdlist: '%s'\n", nd->cmdlist);
}

static void destroy_new_dialog (GtkWidget *w, new_dialog *nd) 
{
    gtk_main_quit();
    free(nd);
}

static char *est_str (int cmdnum)
{
    switch (cmdnum) {
    case OLS:
	return "OLS";
    case HCCM:
	return "HCCM";
    case HSK:
	return "Heteroskedasticity corrected";
    case CORC:
	return "Cochrane-Orcutt";
    case HILU:
	return "Hildreth-Lu";
    case LOGIT:
	return "Logit";
    case PROBIT:
	return "Probit";
    case POOLED:
	return "Pooled OLS";
    case WLS:
	return "Weighted least squares";
    case TSLS:
	return "Two-stage least squares";
    case AR:
	return "Autoregressive";
    case VAR:
	return "VAR";
    default:
	return "";
    }
}

static char *extra_string (int cmdnum)
{
    switch (cmdnum) {
    case WLS:
	return "Weight variable";
    case TSLS:
	return "Instruments";
    case AR:
	return "List of AR lags";
    case VAR:
	return "Lag order";
    default:
	return "";
    }
}

static void 
dialog_select_row (GtkCList *clist, gint row, gint column, 
		   GdkEventButton *event, new_dialog *nd) 
{
    
    if (event != NULL && event->type == GDK_2BUTTON_PRESS) {
	set_dependent_var (row, nd);
    }
}

static gint
dialog_right_click (GtkWidget *widget, GdkEventButton *event, 
		    new_dialog *nd)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(nd->varlist);
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 
    if (mods & GDK_BUTTON3_MASK) {
	add_independent_callback (NULL, nd);
    }
    return TRUE;
}


void new_edit_dialog (const char *title, const char *oktxt, 
		      void (*okfunc)(), guint cmdcode, int *list) 
{
    GtkWidget *right_vbox, *tempwid;
    GtkWidget *big_hbox, *depvar_hbox, *indepvar_hbox;
    GtkWidget *button_vbox, *scroller;
    new_dialog *nd;
    char modelstr[48];
    int i;
    static int default_var;

    nd = mymalloc(sizeof *nd);
    if (nd == NULL) return;

    nd->default_var = &default_var;
    nd->code = cmdcode;

    nd->dlg = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(nd->dlg), title);

    gtk_signal_connect (GTK_OBJECT (nd->dlg), "destroy", 
			GTK_SIGNAL_FUNC (destroy_new_dialog), 
			nd);    

    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(nd->dlg)->vbox), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(nd->dlg)->vbox), 5);

    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(nd->dlg)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 5);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), TRUE);

    gtk_window_set_position(GTK_WINDOW(nd->dlg), GTK_WIN_POS_MOUSE);
    /* gtk_window_set_default_size(GTK_WINDOW(nd->dlg), 363, 380); */

    sprintf(modelstr, "%s model", est_str(cmdcode));
    tempwid = gtk_label_new(modelstr);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->vbox), 
		       tempwid, FALSE, FALSE, 5);
    gtk_widget_show(tempwid);

    /* the following encloses LHS varlist, depvar and indepvar stuff */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* LHS: list of vars to choose from */
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    nd->varlist = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(nd->varlist));
    for (i=0; i<datainfo->v; i++) {
	gchar *row[2];
	gchar id[5];

        if (hidden_var(i, datainfo)) continue;
	sprintf(id, "%d", i);
	row[0] = id;
	row[1] = datainfo->varname[i];
        gtk_clist_append(GTK_CLIST(nd->varlist), row);
    }
    gtk_clist_set_column_width (GTK_CLIST(nd->varlist), 1, 80);
    /* gtk_clist_set_column_visibility (GTK_CLIST(nd->varlist), 0, FALSE); */
    gtk_clist_set_selection_mode (GTK_CLIST(nd->varlist),
				  GTK_SELECTION_EXTENDED);
    gtk_signal_connect_after (GTK_OBJECT (nd->varlist), "select_row", 
                              GTK_SIGNAL_FUNC(dialog_select_row), nd);
    gtk_signal_connect(GTK_OBJECT(nd->varlist),
		       "button_press_event",
		       (GtkSignalFunc) dialog_right_click, nd);
    gtk_widget_show(nd->varlist); 
    gtk_container_add(GTK_CONTAINER(scroller), nd->varlist);
    gtk_widget_show(scroller);
    gtk_box_pack_start(GTK_BOX(big_hbox), scroller, TRUE, TRUE, 0);

    /* RHS: vertical holder for depvar (top) and indepvars (bottom) */
    right_vbox = gtk_vbox_new(FALSE, 5);

    tempwid = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
    gtk_widget_show(tempwid);

    /* top right: dependent variable */
    tempwid = gtk_label_new("Dependent variable");
    gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
    gtk_widget_show(tempwid);

    depvar_hbox = gtk_hbox_new(FALSE, 5); 

    tempwid = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(depvar_hbox), tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(select_dependent_callback), nd);
    gtk_widget_show(tempwid); 

    nd->depvar = gtk_entry_new_with_max_length(8);
    if (default_var) {
	gtk_entry_set_text(GTK_ENTRY(nd->depvar), 
			   datainfo->varname[default_var]);
	gtk_object_set_user_data(GTK_OBJECT(nd->depvar),
				 GINT_TO_POINTER(default_var));
    }
    gtk_box_pack_start(GTK_BOX(depvar_hbox), nd->depvar, FALSE, FALSE, 0);
    gtk_widget_show(nd->depvar); 

    gtk_box_pack_start(GTK_BOX(right_vbox), depvar_hbox, FALSE, FALSE, 0);
    gtk_widget_show(depvar_hbox); 

    nd->default_check = gtk_check_button_new_with_label("Set as default");
    gtk_box_pack_start(GTK_BOX(right_vbox), nd->default_check, FALSE, TRUE, 0);
    gtk_widget_show(nd->default_check); 

    /* separator */
    tempwid = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
    gtk_widget_show(tempwid);

    /* middle right: used for some estimators */
    if (cmdcode == WLS || cmdcode == AR || cmdcode == TSLS || cmdcode == VAR) {

	tempwid = gtk_label_new(extra_string(cmdcode));
	gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
	gtk_widget_show(tempwid);	

	if (cmdcode == WLS) {
	    GtkWidget *midhbox;

	    midhbox = gtk_hbox_new(FALSE, 5);
	    tempwid = gtk_button_new_with_label (_("Choose ->"));
	    gtk_box_pack_start(GTK_BOX(midhbox), tempwid, TRUE, TRUE, 0);
	    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
			    GTK_SIGNAL_FUNC(select_weight_callback), nd);
	    gtk_widget_show(tempwid); 

	    nd->extra_entry = gtk_entry_new_with_max_length(8);
	    gtk_box_pack_start(GTK_BOX(midhbox), nd->extra_entry, FALSE, TRUE, 0);
	    gtk_widget_show(nd->extra_entry); 

	    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, FALSE, TRUE, 0);
	    gtk_widget_show(midhbox); 
	} 
	else if (cmdcode == AR) {
	    nd->extra_entry = gtk_entry_new_with_max_length(8);
	    gtk_box_pack_start(GTK_BOX(right_vbox), nd->extra_entry, 
			       FALSE, TRUE, 0);
	    gtk_widget_show(nd->extra_entry); 
	}
	/* need to deal with other cases: TSLS, VAR */

	/* separator */
	tempwid = gtk_hseparator_new();
	gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
	gtk_widget_show(tempwid);
    }
    
    /* lower right: independent variables */
    tempwid = gtk_label_new("Independent variables");
    gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
    gtk_widget_show(tempwid);

    indepvar_hbox = gtk_hbox_new(FALSE, 5);

    /* push/pull buttons first, in their own little vbox */
    button_vbox = gtk_vbox_new(TRUE, 5);

    tempwid = gtk_button_new_with_label (_("Add ->"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tempwid, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(add_independent_callback), nd);
    gtk_widget_show(tempwid);
    
    tempwid = gtk_button_new_with_label (_("<- Remove"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tempwid, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(remove_independent_callback), nd);
    gtk_widget_show(tempwid);

    gtk_box_pack_start(GTK_BOX(indepvar_hbox), button_vbox, TRUE, TRUE, 0);
    gtk_widget_show(button_vbox);

    /* then the listing */
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    nd->indepvars = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(nd->indepvars));
    for (i=0; i<1; i++) {
	gchar *row[2];
	gchar id[4];

	sprintf(id, "%d", i);
	row[0] = id;
	row[1] = datainfo->varname[i];
        gtk_clist_append(GTK_CLIST(nd->indepvars), row);
    }
    gtk_clist_set_column_width (GTK_CLIST(nd->indepvars), 1, 80);
    gtk_widget_set_usize (nd->indepvars, 80, 120);
    /* gtk_clist_set_column_visibility (GTK_CLIST(nd->indepvars), 0, FALSE); */
    gtk_widget_show(nd->indepvars); 
    gtk_container_add(GTK_CONTAINER(scroller), nd->indepvars);

    gtk_widget_show(scroller);
    gtk_box_pack_start(GTK_BOX(indepvar_hbox), scroller, TRUE, TRUE, 0);

    /* pack the lower right stuff into the RHS vbox */
    gtk_box_pack_start(GTK_BOX(right_vbox), indepvar_hbox, TRUE, TRUE, 0);
    gtk_widget_show(indepvar_hbox);

    /* pack the whole RHS to the right of the LHS varlist */
    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, FALSE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->vbox), 
		       big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    /* buttons: "OK", Clear, Cancel, Help */
    tempwid = gtk_button_new_with_label (oktxt);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(construct_cmdlist), nd);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(okfunc), nd);
    gtk_signal_connect_object(GTK_OBJECT (tempwid), "clicked", 
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      GTK_OBJECT(nd->dlg));
    gtk_widget_show(tempwid);
    gtk_widget_grab_default(tempwid);

    tempwid = gtk_button_new_with_label(_("Clear"));
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(clear_vars), nd);
    gtk_widget_show(tempwid);

    tempwid = gtk_button_new_with_label(_("Cancel"));
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
                       GTK_SIGNAL_FUNC(delete_widget), nd->dlg);
    gtk_widget_show(tempwid);

    tempwid = gtk_button_new_with_label(_("Help"));
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT (tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(context_help), 
		       GINT_TO_POINTER(cmdcode));
    gtk_widget_show(tempwid);

    gtk_widget_show(nd->dlg);
    gtk_main();
}
