#include "gretl.h"

typedef struct {
    GtkWidget *dlg;
    GtkWidget *varlist;
    GtkWidget *depvar;
    GtkWidget *indepvars;
    char cmdlist[MAXLEN];
} new_dialog;

static void set_dependent_var (gint i, new_dialog *nd)
{
    gchar *vname;

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &vname);
    gtk_entry_set_text(GTK_ENTRY(nd->depvar), vname);
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
    gchar *row[1];

    gtk_clist_get_text(GTK_CLIST(nd->varlist), i, 0, &row[0]);
    gtk_clist_append(GTK_CLIST(nd->indepvars), row);
}

static void add_independent_callback (GtkWidget *w, new_dialog *nd)
{
    GList *mylist = GTK_CLIST(nd->varlist)->selection;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) add_independent_var, nd);
    }
}

static void remove_independent_var (gint i, new_dialog *nd)
{
    gtk_clist_remove(GTK_CLIST(nd->indepvars), i);
}

static void remove_independent_callback (GtkWidget *w, new_dialog *nd)
{
    GList *mylist = GTK_CLIST(nd->indepvars)->selection;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) remove_independent_var, nd);
    }
}

static void clear_vars (GtkWidget *w, new_dialog *nd)
{
    gtk_entry_set_text(GTK_ENTRY(nd->depvar), "");    
    gtk_clist_clear(GTK_CLIST(nd->indepvars));
}

static void construct_cmdlist (GtkWidget *w, new_dialog *nd)
{
    nd->cmdlist[0] = 0;

    strcat(nd->cmdlist, gtk_entry_get_text(GTK_ENTRY(nd->depvar)));
    /* now concat all the contents of nd->indepvars */
}

void new_edit_dialog (const char *title, const char *oktxt, 
		      void (*okfunc)(), guint cmdcode) 
{
    GtkWidget *right_vbox, *tempwid;
    GtkWidget *big_hbox, *depvar_hbox, *indepvar_hbox;
    GtkWidget *button_vbox, *scroller;
    new_dialog *nd;
    int i;

    nd = mymalloc(sizeof *nd);
    if (nd == NULL) return;

    nd->dlg = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(nd->dlg), title);
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(nd->dlg)->vbox), 10);
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(nd->dlg)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(nd->dlg)->vbox), 5);

    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(nd->dlg), GTK_WIN_POS_MOUSE);
    /* gtk_window_set_default_size(GTK_WINDOW(nd->dlg), 363, 380); */

    tempwid = gtk_label_new("Select model specification");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->vbox), tempwid, FALSE, FALSE, 5);
    gtk_widget_show(tempwid);

    /* the following encloses LHS varlist, depvar and indepvar stuff */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* LHS: list of vars to choose from */
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

    nd->varlist = gtk_clist_new(1);
    gtk_clist_clear(GTK_CLIST(nd->varlist));
    for (i=0; i<datainfo->v; i++) {
	gchar *row[1];
        if (hidden_var(i, datainfo)) continue;
	row[0] = datainfo->varname[i];
        gtk_clist_append(GTK_CLIST(nd->varlist), row);
    }
    gtk_clist_set_column_width (GTK_CLIST(nd->varlist), 0, 80);
    gtk_widget_show(nd->varlist);  

    gtk_container_add(GTK_CONTAINER(scroller), nd->varlist);
    gtk_box_pack_start(GTK_BOX(big_hbox), scroller, TRUE, TRUE, 0);
    gtk_widget_show(scroller);

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

    nd->depvar = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(nd->depvar), 8);
    gtk_box_pack_start(GTK_BOX(depvar_hbox), nd->depvar, FALSE, FALSE, 0);
    gtk_widget_show(nd->depvar); 
    gtk_box_pack_start(GTK_BOX(right_vbox), depvar_hbox, FALSE, FALSE, 0);
    gtk_widget_show(depvar_hbox); 

    tempwid = gtk_check_button_new_with_label("Set as default");
    gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
    gtk_widget_show(tempwid); 

    /* separator: depvar versus indepvars */
    tempwid = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tempwid, FALSE, TRUE, 0);
    gtk_widget_show(tempwid);
    
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

    nd->indepvars = gtk_clist_new(1);
    gtk_clist_clear(GTK_CLIST(nd->indepvars));
    gtk_clist_set_column_width (GTK_CLIST(nd->indepvars), 0, 80);
    gtk_widget_show(nd->indepvars); 

    gtk_container_add(GTK_CONTAINER(scroller), nd->indepvars);
    gtk_widget_show(scroller);    
    gtk_box_pack_start(GTK_BOX(indepvar_hbox), scroller, TRUE, TRUE, 0);
    gtk_widget_show(indepvar_hbox);

    /* pack the lower right stuff into the RHS vbox */
    gtk_box_pack_start(GTK_BOX(right_vbox), indepvar_hbox, TRUE, TRUE, 0);
    gtk_widget_show(indepvar_hbox);

    /* pack the whole RHS to the right of the LHS varlist */
    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, FALSE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->vbox), big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    /* buttons: "OK", Clear, Cancel, Help */
    tempwid = gtk_button_new_with_label (oktxt);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(okfunc), nd);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label(_("Clear"));
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(clear_vars), nd);
    gtk_widget_show(tempwid);

    tempwid = gtk_button_new_with_label(_("Cancel"));
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
                       GTK_SIGNAL_FUNC(delete_widget), nd->dlg);
    gtk_widget_show(tempwid);

    tempwid = gtk_button_new_with_label(_("Help"));
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nd->dlg)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(context_help), 
			GINT_TO_POINTER(cmdcode));
    gtk_widget_show (tempwid);

    gtk_widget_show (nd->dlg);
}
