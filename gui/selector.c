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
#include "selector.h"

extern GtkWidget *open_dialog; /* dialogs.c */

static int default_var;
static int *xlist;
static int *instlist;
static GtkWidget *scatters_label;
static GtkWidget *scatters_menu;
static GtkWidget *x_axis_item;

void clear_selector (void)
{
    default_var = 0;
    free(xlist);
    xlist = NULL;
    free(instlist);
    instlist = NULL;
}

static gint list_sorter (gconstpointer a, gconstpointer b)
{
    return GPOINTER_TO_INT(b) - GPOINTER_TO_INT(a);
}

static void set_extra_var (gint i, selector *sr)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(sr->extra), vname);
    gtk_object_set_user_data(GTK_OBJECT(sr->extra),
			     GINT_TO_POINTER(atoi(vnum)));
}

static void select_extra_var_callback (GtkWidget *w, selector *sr)
{
    GList *mylist;

    if (!GTK_IS_CLIST(sr->varlist)) return;

    mylist = GTK_CLIST(sr->varlist)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) set_extra_var, sr);
    }
}

static void set_dummy (gint i, selector *sr)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(sr->rightvars), vname);
    gtk_object_set_user_data(GTK_OBJECT(sr->rightvars),
			     GINT_TO_POINTER(atoi(vnum)));
}

static void select_factor_callback (GtkWidget *w, selector *sr)
{
    GList *mylist = GTK_CLIST(sr->varlist)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) set_dummy, sr);
    }
}

static void set_dependent_var (gint i, selector *sr)
{
    gchar *vnum, *vname;

    if (sr->depvar == NULL) return;

    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vname);
    gtk_object_set_user_data(GTK_OBJECT(sr->depvar),
			     GINT_TO_POINTER(atoi(vnum)));
}

static void select_dependent_callback (GtkWidget *w, selector *sr)
{
    GList *mylist;

    if (!GTK_IS_CLIST(sr->varlist)) return;

    mylist = GTK_CLIST(sr->varlist)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) set_dependent_var, sr);
    }
}

static void add_instrument (gint i, selector *sr)
{
    gchar *row[2];
    gint j, rows = GTK_CLIST(sr->extra)->rows;
    gint already_there = 0;

    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 0, &row[0]);
    for (j=0; j<rows; j++) {
	gchar *test;

	gtk_clist_get_text(GTK_CLIST(sr->extra), j, 0, &test);
	if (!strcmp(test, row[0])) {
	    already_there = 1; 
	    break;
	}
    }
    if (!already_there) {
	gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 1, &row[1]);
	gtk_clist_append(GTK_CLIST(sr->extra), row);
    }
}

static void add_instrument_callback (GtkWidget *w, selector *sr)
{
    GList *mylist = GTK_CLIST(sr->varlist)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) add_instrument, sr);
}

static void add_var_on_right (gint i, selector *sr)
{
    gchar *row[2];
    gint j, rows;
    gint already_there = 0;

    if (!GTK_IS_CLIST(sr->rightvars)) return;

    rows = GTK_CLIST(sr->rightvars)->rows;

    gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 0, &row[0]);
    for (j=0; j<rows; j++) {
	gchar *test;

	gtk_clist_get_text(GTK_CLIST(sr->rightvars), j, 0, &test);
	if (!strcmp(test, row[0])) {
	    already_there = 1; 
	    break;
	}
    }
    if (!already_there) {
	gtk_clist_get_text(GTK_CLIST(sr->varlist), i, 1, &row[1]);
	gtk_clist_append(GTK_CLIST(sr->rightvars), row);
    }
}

static void add_to_right_callback (GtkWidget *w, selector *sr)
{
    GList *mylist;

    if (!GTK_IS_CLIST(sr->varlist) ||
	!GTK_IS_CLIST(sr->rightvars)) return;

    mylist = GTK_CLIST(sr->varlist)->selection;

    if (mylist != NULL) 
	g_list_foreach(mylist, (GFunc) add_var_on_right, sr);
}

static void remove_right_var (gint i, selector *sr)
{
    gtk_clist_remove(GTK_CLIST(sr->rightvars), i);
}

static void remove_from_right_callback (GtkWidget *w, selector *sr)
{
    GList *mylist = g_list_copy(GTK_CLIST(sr->rightvars)->selection);
    mylist = g_list_sort(mylist, list_sorter);

    g_list_foreach(mylist, (GFunc) remove_right_var, sr);
}

static void remove_instrument (gint i, selector *sr)
{
    gtk_clist_remove(GTK_CLIST(sr->extra), i);
}

static void remove_instrument_callback (GtkWidget *w, selector *sr)
{
    GList *mylist = g_list_copy(GTK_CLIST(sr->extra)->selection);
    mylist = g_list_sort(mylist, list_sorter);

    g_list_foreach(mylist, (GFunc) remove_instrument, sr);
}

static void clear_vars (GtkWidget *w, selector *sr)
{
    gchar *row[2];

    gtk_clist_unselect_all(GTK_CLIST(sr->varlist));
    if (sr->depvar != NULL) 
	gtk_entry_set_text(GTK_ENTRY(sr->depvar), "");
    if (sr->code == GR_DUMMY)
	gtk_entry_set_text(GTK_ENTRY(sr->rightvars), "");
    else
	gtk_clist_clear(GTK_CLIST(sr->rightvars));
    if (MODEL_CODE(sr->code)) {
	row[0] = "0";
	row[1] = "const";
	gtk_clist_append(GTK_CLIST(sr->rightvars), row);
    }
}

static void topslot_empty (int code)
{
    switch (code) {
    case GR_XY:
    case GR_IMP:
	errbox(_("You must select an X-axis variable"));
	break;
    case SCATTERS:
	errbox(_("You must select a Y-axis variable"));
	break;
    default:
	errbox(_("You must select a dependent variable"));
    }
}

static void reverse_list (char *list)
{
    char *tmp, *p;
    char istr[8];

    p = strchr(list, ';');
    if (p == NULL) return;

    tmp = malloc(strlen(list) + 4);
    if (tmp == NULL) return;

    sscanf(list, "%7s", istr);

    strcpy(tmp, p + 2);
    strcat(tmp, " ; ");
    strcat(tmp, istr);

    strcpy(list, tmp);
    
    free(tmp);
}

static void construct_cmdlist (GtkWidget *w, selector *sr)
{
    gint i = 0, rows = 0;
    gchar numstr[6], grvar[6];
    int err = 0;

    sr->cmdlist = mymalloc(MAXLEN);
    if (sr->cmdlist == NULL) return;
    sr->cmdlist[0] = 0;

    if (sr->code != GR_DUMMY)
	rows = GTK_CLIST(sr->rightvars)->rows;

    /* first deal with content of "extra" widget */
    if (sr->code == WLS) {
	gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->extra));

	if (str == NULL || !strlen(str)) {
	    errbox(_("You must select a weight variable"));
	    err = 1;
	} else {
	    i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(sr->extra)));
	    sprintf(numstr, "%d ", i);
	    strcat(sr->cmdlist, numstr);
	}
    }
    else if (sr->code == AR) {
	gchar *lags;

	lags = gtk_entry_get_text(GTK_ENTRY(sr->extra));
	if (!strlen(lags)) {
	    errbox(_("You must specify a list of lags"));
	    err = 1;
	} else {
	    strcat(sr->cmdlist, lags);
	    strcat(sr->cmdlist, " ; ");
	}
    }
    else if (sr->code == VAR || sr->code == COINT || sr->code == COINT2) {
	GtkAdjustment *adj = 
	    gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(sr->extra));

	i = (gint) adj->value;
	sprintf(numstr, "%d ", i);
	strcat(sr->cmdlist, numstr);
    }
    else if (sr->code == GR_DUMMY) {
	gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->extra));

	if (str == NULL || !strlen(str)) {
	    errbox(_("You must select a Y-axis variable"));
	    err = 1;
	} else {
	    i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(sr->extra)));
	    sprintf(numstr, "%d ", i);
	    strcat(sr->cmdlist, numstr);
	}
    }

    /* next deal with the "depvar" widget */
    if (!err && sr->depvar != NULL) {
	gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->depvar));

	if (str == NULL || !strlen(str)) {
	    topslot_empty(sr->code);
	    err = 1;
	} else {
	    i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(sr->depvar)));
	    if (sr->code == GR_XY || sr->code == GR_IMP)
		sprintf(grvar, " %d", i);
	    else {
		sprintf(numstr, "%d", i);
		strcat(sr->cmdlist, numstr);
	    }
	}
    }

    /* bail out if things have gone wrong already */
    if (err) {
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "clicked");
	return;
    }

    if (sr->default_check != NULL && 
	GTK_TOGGLE_BUTTON(sr->default_check)->active) 
	default_var = i;

    if (sr->code == SCATTERS) strcat(sr->cmdlist, ";");

    if (sr->code == GR_DUMMY) { /* special case */
	gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->rightvars));

	if (str == NULL || !strlen(str)) {
	    errbox(_("You must select a factor variable"));
	    err = 1;
	} else {
	    i = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(sr->rightvars)));
	    sprintf(numstr, " %d", i);
	    strcat(sr->cmdlist, numstr);
	}
	if (err) gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "clicked");
	return;
    }

    if (MODEL_CODE(sr->code)) {
	if (rows > 0) { 
	    xlist = realloc(xlist, (rows + 1) * sizeof(int));
	    if (xlist != NULL) xlist[0] = rows;
	}
    }
    for (i=0; i<rows; i++) {
	gchar *rvar;

	gtk_clist_get_text(GTK_CLIST(sr->rightvars), i, 0, &rvar);
	strcat(sr->cmdlist, " ");
	strcat(sr->cmdlist, rvar);
	if (MODEL_CODE(sr->code) && xlist != NULL) 
	    xlist[i+1] = atoi(rvar);
    }

    if (sr->code == TSLS) {
	rows = GTK_CLIST(sr->extra)->rows;
	if (rows > 0) {
	    instlist = realloc(instlist, (rows + 1) * sizeof(int));
	    if (instlist != NULL) instlist[0] = rows;
	    strcat(sr->cmdlist, " ;");
	    for (i=0; i<rows; i++) {
		gchar *inst;

		gtk_clist_get_text(GTK_CLIST(sr->extra), i, 0, &inst);
		strcat(sr->cmdlist, " ");
		strcat(sr->cmdlist, inst);
		if (instlist != NULL) instlist[i+1] = atoi(inst);
	    }
	} else {
	    errbox(_("You must specify a set of instrumental variables"));
	    err = 1;
	}
    }

    if (sr->code == GR_XY || sr->code == GR_IMP)
	strcat(sr->cmdlist, grvar);

    if (sr->code == SCATTERS) {
	GtkWidget *active_item = GTK_OPTION_MENU(scatters_menu)->menu_item;

	if (active_item == x_axis_item) {
	    reverse_list(sr->cmdlist);
	}
    }

    if (err) 
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "clicked");
}

static void destroy_selector (GtkWidget *w, selector *sr) 
{
    gtk_main_quit();
    free(sr->cmdlist);
    free(sr);
    open_dialog = NULL;
}

static char *est_str (int cmdnum)
{
    switch (cmdnum) {
    case OLS:
	return N_("OLS");
    case HCCM:
	return N_("HCCM");
    case HSK:
	return N_("Heteroskedasticity corrected");
    case CORC:
	return N_("Cochrane-Orcutt");
    case HILU:
	return N_("Hildreth-Lu");
    case LOGIT:
	return N_("Logit");
    case PROBIT:
	return N_("Probit");
    case POOLED:
	return N_("Pooled OLS");
    case WLS:
	return N_("Weighted least squares");
    case TSLS:
	return N_("Two-stage least squares");
    case AR:
	return N_("Autoregressive");
    case VAR:
	return N_("VAR");
    case LAD:
	return N_("LAD");
    case COINT:
    case COINT2:
	return N_("Cointegration");
#ifdef ENABLE_GMP
    case MPOLS:
	return N_("High precision OLS");
#endif
    default:
	return "";
    }
}

static char *extra_string (int cmdnum)
{
    switch (cmdnum) {
    case WLS:
	return N_("Weight variable");
    case TSLS:
	return N_("Instruments");
    case AR:
	return N_("List of AR lags");
    case GR_DUMMY:
	return N_("Y-axis variable");
    default:
	return NULL;
    }
}

static void 
dblclick_dialog_row (GtkCList *clist, gint row, gint column, 
		     GdkEventButton *event, selector *sr) 
{
    if (event != NULL && event->type == GDK_2BUTTON_PRESS) { 
	set_dependent_var (row, sr);
	if (sr->default_check != NULL) 
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(sr->default_check),
					  TRUE);
    }
}

static gint
remove_right_click (GtkWidget *widget, GdkEventButton *event, 
		    selector *sr)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(sr->rightvars);
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 
    if (mods & GDK_BUTTON3_MASK) 
	remove_from_right_callback (NULL, sr);
    return TRUE;
}

static gint
dialog_right_click (GtkWidget *widget, GdkEventButton *event, 
		    selector *sr)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(sr->varlist);
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 
    if (mods & GDK_BUTTON3_MASK)  
	add_to_right_callback (NULL, sr);
    return TRUE;
}

static gint flip_scatters_axis (GtkMenuItem *m, GtkOptionMenu *popdown)
{
    GtkWidget *active_item = popdown->menu_item;

    if (active_item == x_axis_item) {
	gtk_label_set_text(GTK_LABEL(scatters_label), _("Y-axis variables"));
    } else {
	gtk_label_set_text(GTK_LABEL(scatters_label), _("X-axis variables"));
    }

    return FALSE;
}

static GtkWidget *
scatters_popdown (void)
{
    GtkWidget *popdown;
    GtkWidget *menu;
    GtkWidget *tmp;
    const char *popstrings[] = {
        N_("Y-axis variable"),
        N_("X-axis variable")
    };
    int i;

    popdown = gtk_option_menu_new();
    menu = gtk_menu_new();

    for (i=0; i<2; i++) {
        tmp = gtk_menu_item_new_with_label(_(popstrings[i]));
	gtk_signal_connect(GTK_OBJECT(tmp), "activate",
			   GTK_SIGNAL_FUNC(flip_scatters_axis), popdown);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), tmp);
	if (i == 1) x_axis_item = tmp;
    }

    gtk_option_menu_set_menu(GTK_OPTION_MENU(popdown), menu);

    scatters_menu = popdown;

    return popdown;
}

static void build_x_axis_section (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *x_hbox;

    if (sr->code == SCATTERS) {
	tmp = scatters_popdown();
	gtk_widget_show_all(tmp);
    } else {
	tmp = gtk_label_new(_("X-axis variable"));
	gtk_widget_show(tmp);
    }

    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);

    x_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(x_hbox), tmp, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
			GTK_SIGNAL_FUNC(select_dependent_callback), sr);
    gtk_widget_show(tmp); 

    sr->depvar = gtk_entry_new_with_max_length(8);

    gtk_box_pack_start(GTK_BOX(x_hbox), sr->depvar, FALSE, FALSE, 0);
    gtk_widget_show(sr->depvar); 

    gtk_box_pack_start(GTK_BOX(right_vbox), x_hbox, FALSE, FALSE, 0);
    gtk_widget_show(x_hbox); 

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);
}

static void build_depvar_section (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *depvar_hbox;

    if (sr->code == VAR)
	tmp = gtk_label_new (_("First dependent variable"));
    else
	tmp = gtk_label_new (_("Dependent variable"));
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);

    depvar_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(depvar_hbox), tmp, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
			GTK_SIGNAL_FUNC(select_dependent_callback), sr);
    gtk_widget_show(tmp); 

    sr->depvar = gtk_entry_new_with_max_length(8);
    if (default_var) {
	gtk_entry_set_text(GTK_ENTRY(sr->depvar), 
			   datainfo->varname[default_var]);
	gtk_object_set_user_data(GTK_OBJECT(sr->depvar),
				 GINT_TO_POINTER(default_var));
    }
    gtk_box_pack_start(GTK_BOX(depvar_hbox), sr->depvar, FALSE, FALSE, 0);
    gtk_widget_show(sr->depvar); 

    gtk_box_pack_start(GTK_BOX(right_vbox), depvar_hbox, FALSE, FALSE, 0);
    gtk_widget_show(depvar_hbox); 

    sr->default_check = gtk_check_button_new_with_label(_("Set as default"));
    gtk_box_pack_start(GTK_BOX(right_vbox), sr->default_check, FALSE, TRUE, 0);
    gtk_widget_show(sr->default_check); 

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);
}

static void lag_order_spin (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *midhbox;
    GtkObject *adj;
    gfloat order = datainfo->pd;

    midhbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("lag order:"));
    adj = gtk_adjustment_new(order, 1, 24, 1, 1, 1);
    sr->extra = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    gtk_box_pack_start (GTK_BOX (midhbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);
    gtk_box_pack_start (GTK_BOX (midhbox), sr->extra, FALSE, FALSE, 5);
    gtk_widget_show(sr->extra);

    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, FALSE, TRUE, 0);
    gtk_widget_show(midhbox); 
}

static void dummy_box (selector *sr, GtkWidget *hbox)
{
    GtkWidget *tmp;

    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
			GTK_SIGNAL_FUNC(select_factor_callback), sr);
    gtk_widget_show(tmp); 

    sr->rightvars = gtk_entry_new_with_max_length(8);
    gtk_box_pack_start(GTK_BOX(hbox), sr->rightvars, FALSE, TRUE, 0);
    gtk_widget_show(sr->rightvars); 
}

static void extra_var_box (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *midhbox;

    midhbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(midhbox), tmp, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
			GTK_SIGNAL_FUNC(select_extra_var_callback), sr);
    gtk_widget_show(tmp); 

    sr->extra = gtk_entry_new_with_max_length(8);
    gtk_box_pack_start(GTK_BOX(midhbox), sr->extra, FALSE, TRUE, 0);
    gtk_widget_show(sr->extra); 

    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, FALSE, TRUE, 0);
    gtk_widget_show(midhbox); 
}

static void tsls_box (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *midhbox, *button_vbox;
    GtkWidget *scroller;

    midhbox = gtk_hbox_new(FALSE, 5);

    button_vbox = gtk_vbox_new(TRUE, 5);

    tmp = gtk_button_new_with_label (_("Add ->"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(add_instrument_callback), sr);
    gtk_widget_show(tmp);
    
    tmp = gtk_button_new_with_label (_("<- Remove"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(remove_instrument_callback), sr);
    gtk_widget_show(tmp);

    gtk_box_pack_start(GTK_BOX(midhbox), button_vbox, TRUE, TRUE, 0);
    gtk_widget_show(button_vbox);

    /* then the listing */
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    sr->extra = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(sr->extra));

    if (instlist != NULL) {
	int i;

	for (i=1; i<=instlist[0]; i++) {
	    gchar *row[2];
	    gchar id[4];

	    sprintf(id, "%d", instlist[i]);
	    row[0] = id;
	    row[1] = datainfo->varname[instlist[i]];
	    gtk_clist_append(GTK_CLIST(sr->extra), row);
	}
    } else {
	gchar *row[2];

	row[0] = "0";
	row[1] = "const";
	gtk_clist_append(GTK_CLIST(sr->extra), row);
    }

    gtk_clist_set_column_width (GTK_CLIST(sr->extra), 1, 80 * gui_scale);
    gtk_widget_set_usize (sr->extra, 80 * gui_scale, 120 * gui_scale);

    gtk_widget_show(sr->extra); 
    gtk_container_add(GTK_CONTAINER(scroller), sr->extra);

    gtk_widget_show(scroller);
    gtk_box_pack_start(GTK_BOX(midhbox), scroller, TRUE, TRUE, 0);

    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, FALSE, TRUE, 0);
    gtk_widget_show(midhbox); 
}

static void build_mid_section (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp;
    const char *str = _(extra_string(sr->code));

    if (str != NULL) {
	tmp = gtk_label_new(str);
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
	gtk_widget_show(tmp);
    }	

    if (sr->code == WLS || sr->code == GR_DUMMY) 
	extra_var_box (sr, right_vbox);
    else if (sr->code == VAR || sr->code == COINT || sr->code == COINT2)
	lag_order_spin (sr, right_vbox);
    else if (sr->code == TSLS)
	tsls_box (sr, right_vbox);
    else if (sr->code == AR) {
	sr->extra = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(right_vbox), sr->extra, 
			   FALSE, TRUE, 0);
	gtk_widget_show(sr->extra); 
    }

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);
}

static int screen_scalar (int i, int c)
{
    if ((MODEL_CODE(c) || GRAPH_CODE(c) || c == LAGS || c == DIFF || c == LDIFF)
	&& datainfo->vector[i] == 0)
	return 1;
    return 0;
}

static void selector_init (selector *sr, guint code, const char *title)
{
    sr->varlist = NULL;
    sr->depvar = NULL;
    sr->rightvars = NULL;
    sr->default_check = NULL;
    sr->extra = NULL;
    sr->cmdlist = NULL;
    sr->data = NULL;

    sr->code = code;
    sr->dlg = gtk_dialog_new();
    open_dialog = sr->dlg;
    gtk_window_set_title(GTK_WINDOW(sr->dlg), title);

    gtk_signal_connect (GTK_OBJECT (sr->dlg), "destroy", 
			GTK_SIGNAL_FUNC (destroy_selector), 
			sr);    

    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(sr->dlg)->vbox), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(sr->dlg)->vbox), 5);

    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(sr->dlg)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(sr->dlg)->action_area), 5);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(sr->dlg)->action_area), TRUE);

    gtk_window_set_position(GTK_WINDOW(sr->dlg), GTK_WIN_POS_MOUSE);
}    

static void 
build_selector_buttons (selector *sr, const char *oktxt, void (*okfunc)())
{
    GtkWidget *tmp;

    tmp = gtk_button_new_with_label(oktxt);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->action_area), 
		       tmp, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked", 
		       GTK_SIGNAL_FUNC(construct_cmdlist), sr);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked", 
		       GTK_SIGNAL_FUNC(okfunc), sr);
    gtk_signal_connect_object(GTK_OBJECT (tmp), "clicked", 
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      GTK_OBJECT(sr->dlg));
    gtk_widget_show(tmp);
    gtk_widget_grab_default(tmp);

    tmp = gtk_button_new_with_label(_("Clear"));
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->action_area), 
		       tmp, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked", 
		       GTK_SIGNAL_FUNC(clear_vars), sr);
    gtk_widget_show(tmp);

    tmp = gtk_button_new_with_label(_("Cancel"));
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->action_area), 
		       tmp, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
                       GTK_SIGNAL_FUNC(delete_widget), sr->dlg);
    gtk_widget_show(tmp);

    if (sr->code != PRINT && !SAVE_DATA_ACTION(sr->code)) {
	tmp = gtk_button_new_with_label(_("Help"));
	GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->action_area), 
			   tmp, TRUE, TRUE, 0);
	gtk_signal_connect(GTK_OBJECT (tmp), "clicked", 
			   GTK_SIGNAL_FUNC(context_help), 
			   GINT_TO_POINTER(sr->code));
	gtk_widget_show(tmp);
    }
}

void selection_dialog (const char *title, const char *oktxt, 
		       void (*okfunc)(), guint cmdcode) 
{
    GtkWidget *right_vbox, *tmp;
    GtkWidget *big_hbox, *indepvar_hbox;
    GtkWidget *button_vbox, *scroller;
    selector *sr;
    char topstr[48];
    int i;

    if (open_dialog != NULL) {
	gdk_window_raise(open_dialog->window);
	return;
    }

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) return;

    selector_init(sr, cmdcode, title);

    if (MODEL_CODE(cmdcode))
	strcpy(topstr, _(est_str(cmdcode)));
    else if (cmdcode == GR_XY)
	strcpy(topstr, _("XY scatterplot"));
    else if (cmdcode == GR_IMP)
	strcpy(topstr, _("plot with impulses"));
    else if (cmdcode == SCATTERS)
	strcpy(topstr, _("multiple scatterplots"));
    else if (cmdcode == GR_DUMMY)
	strcpy(topstr, _("factorized plot"));
    else
	strcpy(topstr, "fixme need string");

    tmp = gtk_label_new(topstr);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->vbox), 
		       tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    /* the following encloses LHS varlist, depvar and indepvar stuff */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* LHS: list of vars to choose from */
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    sr->varlist = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(sr->varlist));
    for (i=0; i<datainfo->v; i++) {
	gchar *row[2];
	gchar id[5];

	if (i == 0 && !MODEL_CODE(cmdcode)) continue;
        if (hidden_var(i, datainfo)) continue;
	if (screen_scalar(i, cmdcode)) continue;
	sprintf(id, "%d", i);
	row[0] = id;
	row[1] = datainfo->varname[i];
        gtk_clist_append(GTK_CLIST(sr->varlist), row);
    }
    gtk_clist_set_column_width (GTK_CLIST(sr->varlist), 1, 80);
    gtk_clist_set_selection_mode (GTK_CLIST(sr->varlist),
				  GTK_SELECTION_EXTENDED);
    gtk_signal_connect_after (GTK_OBJECT (sr->varlist), "select_row", 
			      GTK_SIGNAL_FUNC(dblclick_dialog_row), sr);
    gtk_signal_connect(GTK_OBJECT(sr->varlist), "button_press_event",
		       (GtkSignalFunc) dialog_right_click, sr);
    gtk_widget_show(sr->varlist); 
    gtk_container_add(GTK_CONTAINER(scroller), sr->varlist);
    gtk_widget_show(scroller);
    gtk_box_pack_start(GTK_BOX(big_hbox), scroller, TRUE, TRUE, 0);

    /* RHS: vertical holder */
    right_vbox = gtk_vbox_new(FALSE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);

    /* for models: top right -> dependent variable */
    if (MODEL_CODE(cmdcode)) 
	build_depvar_section(sr, right_vbox);
    /* graphs: top right -> x-axis variable */
    else if (cmdcode == GR_XY || cmdcode == GR_IMP || cmdcode == GR_DUMMY
	     || cmdcode == SCATTERS)
	build_x_axis_section(sr, right_vbox);

    /* middle right: used for some estimators and factored plot */
    if (cmdcode == WLS || cmdcode == AR || cmdcode == TSLS || 
	cmdcode == VAR || cmdcode == COINT || cmdcode == COINT2 ||
	cmdcode == GR_DUMMY) 
	build_mid_section(sr, right_vbox);
    
    /* lower right: selected (independent) variables */
    if (MODEL_CODE(cmdcode)) {
	tmp = gtk_label_new(_("Independent variables"));
    } else if (cmdcode == GR_XY || cmdcode == GR_IMP) {
	tmp = gtk_label_new(_("Y-axis variables"));
    } else if (cmdcode == SCATTERS) {
	scatters_label = tmp = gtk_label_new(_("X-axis variables"));
    } else if (cmdcode == GR_DUMMY) {
	tmp = gtk_label_new(_("Factor (dummy)"));
    }
    
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, TRUE, 0);
    gtk_widget_show(tmp);

    indepvar_hbox = gtk_hbox_new(FALSE, 5);

    if (cmdcode == GR_DUMMY) {
	dummy_box(sr, indepvar_hbox);
    } else { /* all other uses: scrollable list of vars */

	/* push/pull buttons first, in their own little vbox */
	button_vbox = gtk_vbox_new(TRUE, 5);

	tmp = gtk_button_new_with_label (_("Add ->"));
	gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
	gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
			    GTK_SIGNAL_FUNC(add_to_right_callback), sr);
	gtk_widget_show(tmp);
    
	tmp = gtk_button_new_with_label (_("<- Remove"));
	gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
	gtk_signal_connect_after (GTK_OBJECT(tmp), "clicked", 
				  GTK_SIGNAL_FUNC(remove_from_right_callback), sr);
	gtk_widget_show(tmp);

	gtk_box_pack_start(GTK_BOX(indepvar_hbox), button_vbox, TRUE, TRUE, 0);
	gtk_widget_show(button_vbox);

	/* then the listing */
	scroller = gtk_scrolled_window_new (NULL, NULL);
	gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
					GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

	sr->rightvars = gtk_clist_new(2);
	gtk_clist_clear(GTK_CLIST(sr->rightvars));

	if (MODEL_CODE(cmdcode) && xlist != NULL) {
	    for (i=1; i<=xlist[0]; i++) {
		gchar *row[2];
		gchar id[4];

		sprintf(id, "%d", xlist[i]);
		row[0] = id;
		row[1] = datainfo->varname[xlist[i]];
		gtk_clist_append(GTK_CLIST(sr->rightvars), row);
	    }
	} else if (MODEL_CODE(cmdcode) && cmdcode != COINT &&
		   cmdcode != COINT2) {
	    gchar *row[2];

	    row[0] = "0";
	    row[1] = "const";
	    gtk_clist_append(GTK_CLIST(sr->rightvars), row);
	}

	gtk_clist_set_column_width (GTK_CLIST(sr->rightvars), 1, 80 * gui_scale);
	gtk_widget_set_usize (sr->rightvars, 80 * gui_scale, 120 * gui_scale);
	gtk_clist_set_selection_mode (GTK_CLIST(sr->rightvars),
				      GTK_SELECTION_EXTENDED);
	gtk_signal_connect(GTK_OBJECT(sr->rightvars), "button_press_event",
			   (GtkSignalFunc) remove_right_click, sr);
	gtk_widget_show(sr->rightvars); 
	gtk_container_add(GTK_CONTAINER(scroller), sr->rightvars);

	gtk_widget_show(scroller);
	gtk_box_pack_start(GTK_BOX(indepvar_hbox), scroller, TRUE, TRUE, 0);
    }

    /* pack the lower right stuff into the RHS vbox */
    gtk_box_pack_start(GTK_BOX(right_vbox), indepvar_hbox, TRUE, TRUE, 0);
    gtk_widget_show(indepvar_hbox);

    /* pack the whole RHS to the right of the LHS varlist */
    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, FALSE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->vbox), 
		       big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    /* buttons: "OK", Clear, Cancel, Help */
    build_selector_buttons(sr, oktxt, okfunc);

    gtk_widget_show(sr->dlg);
    gtk_main();
}

static char *get_topstr (int cmdnum)
{
    switch (cmdnum) {    
    case LOGS:
	return N_("Select variables for logging");
    case LAGS:
	return N_("Select variables for lagging");
    case SQUARE:
	return N_("Select variables to square");
    case DIFF:
	return N_("Select variables to difference");
    case LDIFF:
	return N_("Select variables to log-difference");
    case ADD:
	return N_("Select variables to add");
    case OMIT:
	return N_("Select variables to omit");
    case COEFFSUM:
	return N_("Select coefficients to sum");
    case PRINT:
	return N_("Select variables to display");
    case GR_PLOT: 
    case GR_BOX: 
    case GR_NBOX:
	return N_("Select variables to plot");
    default:
	return "";
    }
}

static void add_omit_list (gpointer p, selector *sr)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int i;

    if (sr->code == OMIT || sr->code == COEFFSUM) {
	for (i=2; i<=pmod->list[0]; i++) {
	    gchar *row[2];
	    gchar id[5];

	    if (pmod->list[i] == 0) continue;
	    sprintf(id, "%d", pmod->list[i]);
	    row[0] = id;
	    row[1] = datainfo->varname[pmod->list[i]];
	    gtk_clist_append(GTK_CLIST(sr->varlist), row);
	} 
    } else {
	for (i=1; i<datainfo->v; i++) {
	    gchar *row[2];
	    gchar id[5];
	    int j, match = 0;

	    for (j=1; j<=pmod->list[0]; j++) {
		if (i == pmod->list[j]) {
		    match = 1;
		    break;
		}
	    }
	    if (match) continue;
	    sprintf(id, "%d", i);
	    row[0] = id;
	    row[1] = datainfo->varname[i];
	    gtk_clist_append(GTK_CLIST(sr->varlist), row);
	}
    }
}

static GtkWidget *selection_top_label (int code)
{
    GtkWidget *label = NULL;
    const char *str = get_topstr(code);

    if (strlen(str)) {
	label = gtk_label_new(_(str));
    } 

    return label;
}

void simple_selection (const char *title, const char *oktxt, 
		       void (*okfunc)(), guint cmdcode,
		       gpointer p) 
{
    GtkWidget *left_vbox, *mid_vbox, *right_vbox, *tmp;
    GtkWidget *top_hbox, *big_hbox, *scroller;
    selector *sr;
    int i;

    if (open_dialog != NULL) {
	gdk_window_raise(open_dialog->window);
	return;
    }

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) return;

    selector_init(sr, cmdcode, title);

    sr->data = p;

    tmp = selection_top_label(cmdcode);
    if (tmp != NULL) {
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->vbox), 
			   tmp, TRUE, TRUE, 0);
	gtk_widget_show(tmp);
    }

    /* for titles */
    top_hbox = gtk_hbox_new(FALSE, 0); 
    gtk_box_set_homogeneous(GTK_BOX(top_hbox), TRUE);

    tmp = gtk_label_new(_("Available vars"));
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, TRUE, TRUE, 5);
    gtk_widget_show(tmp);

    tmp = gtk_label_new(" ");
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, TRUE, TRUE, 5);
    gtk_widget_show(tmp);

    tmp = gtk_label_new(_("Selected vars"));
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, TRUE, TRUE, 5);
    gtk_widget_show(tmp);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->vbox), top_hbox, 
		       FALSE, FALSE, 5);
    gtk_widget_show(top_hbox);

    /* the following encloses 3 vboxes */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* holds available var list */
    left_vbox = gtk_vbox_new(FALSE, 5);

    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    sr->varlist = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(sr->varlist));
    if (cmdcode == OMIT || cmdcode == ADD || cmdcode == COEFFSUM) {
        add_omit_list(p, sr);
    } else {
	for (i=1; i<datainfo->v; i++) {
	    gchar *row[2];
	    gchar id[5];

	    if (hidden_var(i, datainfo)) continue;
	    if (screen_scalar(i, cmdcode)) continue;
	    sprintf(id, "%d", i);
	    row[0] = id;
	    row[1] = datainfo->varname[i];
	    gtk_clist_append(GTK_CLIST(sr->varlist), row);
	}
    }

    gtk_clist_set_column_width (GTK_CLIST(sr->varlist), 1, 80 * gui_scale);
    gtk_widget_set_usize (sr->varlist, 80 * gui_scale, 120 * gui_scale);
    gtk_clist_set_selection_mode (GTK_CLIST(sr->varlist),
				  GTK_SELECTION_EXTENDED);
    gtk_signal_connect(GTK_OBJECT(sr->varlist), "button_press_event",
		       (GtkSignalFunc) dialog_right_click, sr);
    gtk_container_add(GTK_CONTAINER(scroller), sr->varlist);
    gtk_widget_show(sr->varlist); 
    gtk_box_pack_start(GTK_BOX(left_vbox), scroller, TRUE, TRUE, 0);
    gtk_widget_show(scroller);

    gtk_box_pack_start(GTK_BOX(big_hbox), left_vbox, TRUE, TRUE, 0);
    gtk_widget_show(left_vbox);
    
    /* middle: vertical holder for push/pull buttons */
    mid_vbox = gtk_vbox_new(FALSE, 5);

    tmp = gtk_button_new_with_label (_("Select ->"));
    gtk_box_pack_start(GTK_BOX(mid_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(add_to_right_callback), sr);
    gtk_widget_show(tmp);
    
    tmp = gtk_button_new_with_label (_("<- Remove"));
    gtk_box_pack_start(GTK_BOX(mid_vbox), tmp, TRUE, FALSE, 0);
    gtk_signal_connect (GTK_OBJECT(tmp), "clicked", 
                        GTK_SIGNAL_FUNC(remove_from_right_callback), sr);
    gtk_widget_show(tmp);

    gtk_box_pack_start(GTK_BOX(big_hbox), mid_vbox, TRUE, TRUE, 0);
    gtk_widget_show(mid_vbox);

    /* RHS: vertical holder for selected vars */
    right_vbox = gtk_vbox_new(FALSE, 5);

    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    sr->rightvars = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(sr->rightvars));

    gtk_clist_set_column_width (GTK_CLIST(sr->rightvars), 1, 80 * gui_scale);
    gtk_widget_set_usize (sr->rightvars, 80 * gui_scale, 120 * gui_scale);
    gtk_clist_set_selection_mode (GTK_CLIST(sr->rightvars),
				  GTK_SELECTION_EXTENDED);
    gtk_container_add(GTK_CONTAINER(scroller), sr->rightvars);
    gtk_widget_show(sr->rightvars); 

    gtk_box_pack_start(GTK_BOX(right_vbox), scroller, TRUE, TRUE, 0);
    gtk_widget_show(scroller);

    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, TRUE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(sr->dlg)->vbox), 
		       big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    /* buttons: "OK", Clear, Cancel, Help */
    build_selector_buttons(sr, oktxt, okfunc);

    gtk_widget_show(sr->dlg);
    gtk_main();
}

static const char *data_save_title (int code)
{
    switch (code) {
    case EXPORT_CSV:
        return _("Save CSV data file");
    case EXPORT_R:
    case EXPORT_R_ALT:
        return _("Save R data file");
    case EXPORT_OCTAVE:
        return _("Save octave data file");
    case COPY_CSV:
        return N_("Select variables to copy");
    default:
        return _("Save data file");
    }
    return "";
}

static void data_save_selection_callback (GtkWidget *w, gpointer p)
{
    selector *sr = (selector *) p;

    if (sr->cmdlist == NULL || *sr->cmdlist == 0) return;

    if (storelist != NULL) {
        free(storelist);
        storelist = NULL;
    }

    storelist = g_strdup(sr->cmdlist);

    if (sr->code != COPY_CSV) {
        file_selector(data_save_title(sr->code), sr->code, NULL);
    }
}

void data_save_selection_wrapper (int file_code)
{
    simple_selection((file_code == COPY_CSV)? 
                     _("Copy data") : _("Save data"), _("Apply"),
                     data_save_selection_callback, file_code, 
                     NULL);
}

struct list_maker {
    char *liststr;
    int n_items;
    size_t len;
    int overflow;
};

static void selection_add_item (gint i, struct list_maker *lmkr)
{
    gchar *varnum = NULL;

    if (lmkr->len > MAXLEN - 12) {
	lmkr->overflow = 1;
	return;
    }

    if (gtk_clist_get_text(GTK_CLIST(mdata->listbox), i, 0, &varnum)) {   
	strcat(lmkr->liststr, " ");
	strcat(lmkr->liststr, varnum);
	lmkr->len += strlen(varnum) + 1;
	lmkr->n_items += 1;
    }
}

char *mdata_selection_to_string (void)
{
    GList *mylist = GTK_CLIST(mdata->listbox)->selection;
    struct list_maker lmkr;    

    lmkr.liststr = mymalloc(MAXLEN);
    if (lmkr.liststr == NULL) return NULL;
    lmkr.liststr[0] = 0;
    lmkr.n_items = lmkr.overflow = 0;
    lmkr.len = 0;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) selection_add_item, 
		       &lmkr);
    }

    if (lmkr.overflow) {
	errbox(_("Too many items were selected"));
	lmkr.liststr[0] = 0;
    }

    return lmkr.liststr;
}
