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

/* selector.c for gretl */

#include "gretl.h"
#include "fancy_dialog.h"

static int default_var;
static int *xlist;
static int *instlist;

static void set_weight_var (gint i, selector *nd)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(nd->extra), vname);
    gtk_object_set_user_data(GTK_OBJECT(nd->extra),
			     GINT_TO_POINTER(atoi(vnum)));
}

static void select_weight_callback (GtkWidget *w, selector *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) set_weight_var, nd);
    }
}

static void set_dependent_var (gint i, selector *nd)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(nd->depvar), vname);
    gtk_object_set_user_data(GTK_OBJECT(nd->depvar),
			     GINT_TO_POINTER(atoi(vnum)));
}

static void select_dependent_callback (GtkWidget *w, selector *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) set_dependent_var, nd);
    }
}

static void add_instrument (gint i, selector *nd)
{
    gchar *row[2];
    gint j, rows = GTK_CLIST(nd->extra)->rows;
    gint already_there = 0;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &row[0]);
    for (j=0; j<rows; j++) {
	gchar *test;

	gtk_clist_get_text(GTK_CLIST(nd->extra), j, 0, &test);
	if (!strcmp(test, row[0])) {
	    already_there = 1; 
	    break;
	}
    }
    if (!already_there) {
	gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 1, &row[1]);
	gtk_clist_append(GTK_CLIST(nd->extra), row);
    }
}

static void add_instrument_callback (GtkWidget *w, selector *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) add_instrument, nd);
}

static void add_independent_var (gint i, selector *nd)
{
    gchar *row[2];
    gint j, rows = GTK_CLIST(nd->rightvars)->rows;
    gint already_there = 0;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &row[0]);
    for (j=0; j<rows; j++) {
	gchar *test;

	gtk_clist_get_text(GTK_CLIST(nd->rightvars), j, 0, &test);
	if (!strcmp(test, row[0])) {
	    already_there = 1; 
	    break;
	}
    }
    if (!already_there) {
	gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 1, &row[1]);
	gtk_clist_append(GTK_CLIST(nd->rightvars), row);
    }
}

static void add_independent_callback (GtkWidget *w, selector *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) add_independent_var, nd);
}

static void remove_independent_var (gint i, selector *nd)
{
    gtk_clist_remove(GTK_CLIST(nd->rightvars), i);
}

static void remove_independent_callback (GtkWidget *w, selector *nd)
{
    GList *mylist = GTK_CLIST(nd->rightvars)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) remove_independent_var, nd);
}

static void remove_instrument (gint i, selector *nd)
{
    gtk_clist_remove(GTK_CLIST(nd->extra), i);
}

static void remove_instrument_callback (GtkWidget *w, selector *nd)
{
    GList *mylist = GTK_CLIST(nd->extra)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) remove_instrument, nd);
}

static void clear_vars (GtkWidget *w, selector *nd)
{
    gchar *row[2];

    gtk_clist_unselect_all(GTK_CLIST(nd->varlist));
    if (nd->depvar != NULL) 
	gtk_entry_set_text(GTK_ENTRY(nd->depvar), ""); 
    gtk_clist_clear(GTK_CLIST(nd->rightvars));
    if (MODEL_CODE(nd->code)) {
	row[0] = "0";
	row[1] = "const";
	gtk_clist_append(GTK_CLIST(nd->rightvars), row);
    }
}

static void construct_cmdlist (GtkWidget *w, selector *nd)
{
    gint i = 0, rows = GTK_CLIST(nd->rightvars)->rows;
    gchar numstr[6];
    int err = 0;

    nd->cmdlist = mymalloc(MAXLEN);
    if (nd->cmdlist == NULL) return;
    nd->cmdlist[0] = 0;

    if (nd->code == WLS) {
	i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(nd->extra)));
	sprintf(numstr, "%d ", i);
	strcat(nd->cmdlist, numstr);
    }
    else if (nd->code == AR) {
	gchar *lags;

	lags = gtk_entry_get_text(GTK_ENTRY(nd->extra));
	strcat(nd->cmdlist, lags);
	strcat(nd->cmdlist, " ; ");
    }
    else if (nd->code == VAR) {
	GtkAdjustment *adj = 
	    gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(nd->extra));

	i = (gint) adj->value;
	sprintf(numstr, "%d ", i);
	strcat(nd->cmdlist, numstr);
    }

    if (nd->depvar != NULL) {
	i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(nd->depvar)));
	sprintf(numstr, "%d", i);
	strcat(nd->cmdlist, numstr);
    }
    
    if (nd->default_check != NULL && GTK_TOGGLE_BUTTON(nd->default_check)->active) 
	default_var = i;

    if (rows > 0) { 
	xlist = realloc(xlist, (rows + 1) * sizeof(int));
	if (xlist != NULL) xlist[0] = rows;
    }

    for (i=0; i<rows; i++) {
	gchar *rvar;

	gtk_clist_get_text(GTK_CLIST(nd->rightvars), i, 0, &rvar);
	strcat(nd->cmdlist, " ");
	strcat(nd->cmdlist, rvar);
	if (xlist != NULL) xlist[i+1] = atoi(rvar);
    }

    if (nd->code == TSLS) {
	rows = GTK_CLIST(nd->extra)->rows;
	if (rows > 0) {
	    instlist = realloc(instlist, (rows + 1) * sizeof(int));
	    if (instlist != NULL) instlist[0] = rows;
	    strcat(nd->cmdlist, " ;");
	    for (i=0; i<rows; i++) {
		gchar *inst;

		gtk_clist_get_text(GTK_CLIST(nd->extra), i, 0, &inst);
		strcat(nd->cmdlist, " ");
		strcat(nd->cmdlist, inst);
		if (instlist != NULL) instlist[i+1] = atoi(inst);
	    }
	} else {
	    errbox("No instrumental variables were specified");
	    err = 1;
	}
    }

    fprintf(stderr, "cmdlist: '%s'\n", nd->cmdlist);
    if (err) 
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "clicked");
}

static void destroy_selector (GtkWidget *w, selector *nd) 
{
    gtk_main_quit();
    fprintf(stderr, "done gtk_main_quit(), now freeing nd\n");
    free(nd->cmdlist);
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

static char *addvar_str (int cmdnum)
{
    switch (cmdnum) {    
    case LOGS:
	return _("for logging");
	break;
    case LAGS:
	return _("for lagging");
    case SQUARE:
	return _("to square");
    case DIFF:
	return _("to difference");
    case LDIFF:
	return _("to log-difference");
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
    default:
	return NULL;
    }
}

static void 
dialog_select_row (GtkCList *clist, gint row, gint column, 
		   GdkEventButton *event, selector *nd) 
{
    
    if (event != NULL && event->type == GDK_2BUTTON_PRESS) {
	set_dependent_var (row, nd);
    }
}

static gint
remove_right_click (GtkWidget *widget, GdkEventButton *event, 
		    selector *nd)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(nd->rightvars);
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 
    if (mods & GDK_BUTTON3_MASK) {
	remove_independent_callback (NULL, nd);
    }
    return TRUE;
}

static gint
dialog_right_click (GtkWidget *widget, GdkEventButton *event, 
		    selector *nd)
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

static void build_depvar_section (selector *nd, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *depvar_hbox;

    if (nd->code == VAR)
	tmp = gtk_label_new("First dependent variable");
    else
	tmp = gtk_label_new("Dependent variable");
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);

    depvar_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(depvar_hbox), tmp, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
			GTK_SIGNAL_FUNC(select_dependent_callback), nd);
    gtk_widget_show(tmp); 

    nd->depvar = gtk_entry_new_with_max_length(9);
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
    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);
}

static void lag_order_spin (selector *nd, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *midhbox;
    GtkObject *adj;
    gfloat order = datainfo->pd;

    midhbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new("lag order:");
    adj = gtk_adjustment_new(order, 1, 24, 1, 1, 1);
    nd->extra = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    gtk_box_pack_start (GTK_BOX (midhbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);
    gtk_box_pack_start (GTK_BOX (midhbox), nd->extra, FALSE, FALSE, 5);
    gtk_widget_show(nd->extra);

    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, FALSE, TRUE, 0);
    gtk_widget_show(midhbox); 
}

static void weight_box (selector *nd, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *midhbox;

    midhbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(midhbox), tmp, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
			GTK_SIGNAL_FUNC(select_weight_callback), nd);
    gtk_widget_show(tmp); 

    nd->extra = gtk_entry_new_with_max_length(8);
    gtk_box_pack_start(GTK_BOX(midhbox), nd->extra, FALSE, TRUE, 0);
    gtk_widget_show(nd->extra); 

    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, FALSE, TRUE, 0);
    gtk_widget_show(midhbox); 
}

static void tsls_box (selector *nd, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *midhbox, *button_vbox;
    GtkWidget *scroller;

    midhbox = gtk_hbox_new(FALSE, 5);

    button_vbox = gtk_vbox_new(TRUE, 5);

    tmp = gtk_button_new_with_label (_("Add ->"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(add_instrument_callback), nd);
    gtk_widget_show(tmp);
    
    tmp = gtk_button_new_with_label (_("<- Remove"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(remove_instrument_callback), nd);
    gtk_widget_show(tmp);

    gtk_box_pack_start(GTK_BOX(midhbox), button_vbox, TRUE, TRUE, 0);
    gtk_widget_show(button_vbox);

    /* then the listing */
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    nd->extra = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(nd->extra));

    if (instlist != NULL) {
	int i;

	for (i=1; i<=instlist[0]; i++) {
	    gchar *row[2];
	    gchar id[4];

	    sprintf(id, "%d", instlist[i]);
	    row[0] = id;
	    row[1] = datainfo->varname[instlist[i]];
	    gtk_clist_append(GTK_CLIST(nd->extra), row);
	}
    } else {
	gchar *row[2];

	row[0] = "0";
	row[1] = "const";
	gtk_clist_append(GTK_CLIST(nd->extra), row);
    }

    gtk_clist_set_column_width (GTK_CLIST(nd->extra), 1, 80);
    gtk_widget_set_usize (nd->extra, 80, 120);

    gtk_widget_show(nd->extra); 
    gtk_container_add(GTK_CONTAINER(scroller), nd->extra);

    gtk_widget_show(scroller);
    gtk_box_pack_start(GTK_BOX(midhbox), scroller, TRUE, TRUE, 0);

    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, FALSE, TRUE, 0);
    gtk_widget_show(midhbox); 
}

static void build_mid_section (selector *nd, GtkWidget *right_vbox)
{
    GtkWidget *tmp;
    char *str = extra_string(nd->code);

    if (str != NULL) {
	tmp = gtk_label_new(str);
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
	gtk_widget_show(tmp);
    }	

    if (nd->code == WLS) 
	weight_box (nd, right_vbox);
    else if (nd->code == VAR)
	lag_order_spin (nd, right_vbox);
    else if (nd->code == TSLS)
	tsls_box (nd, right_vbox);
    else if (nd->code == AR) {
	nd->extra = gtk_entry_new_with_max_length(8);
	gtk_box_pack_start(GTK_BOX(right_vbox), nd->extra, 
			   FALSE, TRUE, 0);
	gtk_widget_show(nd->extra); 
    }

    /* separator */
    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);
}

static void selector_init (selector *nd)
{
    nd->dlg = NULL;
    nd->varlist = NULL;
    nd->depvar = NULL;
    nd->rightvars = NULL;
    nd->default_check = NULL;
    nd->extra = NULL;
    nd->code = 0;
    nd->cmdlist = NULL;
}

void selection_dialog (const char *title, const char *oktxt, 
		       void (*okfunc)(), guint cmdcode) 
{
    GtkWidget *right_vbox, *tmp;
    GtkWidget *big_hbox, *indepvar_hbox;
    GtkWidget *button_vbox, *scroller;
    selector *nd;
    char topstr[48];
    int i;

    nd = mymalloc(sizeof *nd);
    if (nd == NULL) return;
    selector_init(nd);

    nd->code = cmdcode;

    nd->dlg = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(nd->dlg), title);

    gtk_signal_connect (GTK_OBJECT (nd->dlg), "destroy", 
			GTK_SIGNAL_FUNC (destroy_selector), 
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

    if (MODEL_CODE(cmdcode))
	sprintf(topstr, "%s model", est_str(cmdcode));
    else if (ADDVAR_CODE(cmdcode))
	sprintf(topstr, "Select variables %s", addvar_str(cmdcode));
    else
	strcpy(topstr, "fixme need string");
    tmp = gtk_label_new(topstr);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->vbox), 
		       tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

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
    if (MODEL_CODE(cmdcode)) {
	gtk_signal_connect_after (GTK_OBJECT (nd->varlist), "select_row", 
				  GTK_SIGNAL_FUNC(dialog_select_row), nd);
    }
    gtk_signal_connect(GTK_OBJECT(nd->varlist), "button_press_event",
		       (GtkSignalFunc) dialog_right_click, nd);
    gtk_widget_show(nd->varlist); 
    gtk_container_add(GTK_CONTAINER(scroller), nd->varlist);
    gtk_widget_show(scroller);
    gtk_box_pack_start(GTK_BOX(big_hbox), scroller, TRUE, TRUE, 0);

    /* RHS: vertical holder for depvar (top) and indepvars (bottom) */
    right_vbox = gtk_vbox_new(FALSE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);

    /* for models: top right -- dependent variable */
    if (MODEL_CODE(cmdcode)) 
	build_depvar_section(nd, right_vbox);

    /* middle right: used for some estimators */
    if (cmdcode == WLS || cmdcode == AR || cmdcode == TSLS || cmdcode == VAR) 
	build_mid_section(nd, right_vbox);
    
    /* lower right: independent variables */
    if (MODEL_CODE(cmdcode))
	tmp = gtk_label_new("Independent variables");
    else if (ADDVAR_CODE(cmdcode))
	tmp = gtk_label_new("Selected variables");
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);

    indepvar_hbox = gtk_hbox_new(FALSE, 5);

    /* push/pull buttons first, in their own little vbox */
    button_vbox = gtk_vbox_new(TRUE, 5);

    tmp = gtk_button_new_with_label (_("Add ->"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(add_independent_callback), nd);
    gtk_widget_show(tmp);
    
    tmp = gtk_button_new_with_label (_("<- Remove"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(remove_independent_callback), nd);
    gtk_widget_show(tmp);

    gtk_box_pack_start(GTK_BOX(indepvar_hbox), button_vbox, TRUE, TRUE, 0);
    gtk_widget_show(button_vbox);

    /* then the listing */
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    nd->rightvars = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(nd->rightvars));

    if (xlist != NULL) {
	for (i=1; i<=xlist[0]; i++) {
	    gchar *row[2];
	    gchar id[4];

	    sprintf(id, "%d", xlist[i]);
	    row[0] = id;
	    row[1] = datainfo->varname[xlist[i]];
	    gtk_clist_append(GTK_CLIST(nd->rightvars), row);
	}
    } else if (MODEL_CODE(cmdcode)) {
	    gchar *row[2];

	    row[0] = "0";
	    row[1] = "const";
	    gtk_clist_append(GTK_CLIST(nd->rightvars), row);
    }

    gtk_clist_set_column_width (GTK_CLIST(nd->rightvars), 1, 80);
    gtk_widget_set_usize (nd->rightvars, 80, 120);
    gtk_signal_connect(GTK_OBJECT(nd->rightvars), "button_press_event",
		       (GtkSignalFunc) remove_right_click, nd);
    /* gtk_clist_set_column_visibility (GTK_CLIST(nd->rightvars), 0, FALSE); */
    gtk_widget_show(nd->rightvars); 
    gtk_container_add(GTK_CONTAINER(scroller), nd->rightvars);

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
    tmp = gtk_button_new_with_label (oktxt);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tmp, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked", 
		       GTK_SIGNAL_FUNC(construct_cmdlist), nd);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked", 
		       GTK_SIGNAL_FUNC(okfunc), nd);
    gtk_signal_connect_object(GTK_OBJECT (tmp), "clicked", 
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      GTK_OBJECT(nd->dlg));
    gtk_widget_show(tmp);
    gtk_widget_grab_default(tmp);

    tmp = gtk_button_new_with_label(_("Clear"));
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tmp, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked", 
		       GTK_SIGNAL_FUNC(clear_vars), nd);
    gtk_widget_show(tmp);

    tmp = gtk_button_new_with_label(_("Cancel"));
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tmp, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
                       GTK_SIGNAL_FUNC(delete_widget), nd->dlg);
    gtk_widget_show(tmp);

    tmp = gtk_button_new_with_label(_("Help"));
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tmp, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT (tmp), "clicked", 
		       GTK_SIGNAL_FUNC(context_help), 
		       GINT_TO_POINTER(cmdcode));
    gtk_widget_show(tmp);

    gtk_widget_show(nd->dlg);
    gtk_main();
}
