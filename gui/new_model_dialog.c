#include "gretl.h"

void new_model_dialog (gpointer p, guint u, GtkWidget *w) 
{
    GtkWidget *dlg, *hbox, *vbox, *tempwid;
    GtkWidget *varbox, *indepvarbox, *box;
    int i;

    dlg = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(dlg), _("model specification"));
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(dlg)->vbox), 10);
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(dlg)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dlg)->vbox), 2);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dlg)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(dlg)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(dlg), GTK_WIN_POS_MOUSE);

    /* guts of the dialog go here */
    /* want var clist on left, depvar box and indepvar clist on right,
       with buttons to push selected vars (on left) into place */

    hbox = gtk_hbox_new(FALSE, 1);

    varbox = gtk_clist_new(1);
    gtk_clist_clear(GTK_CLIST(varbox));
    for (i=0; i<datainfo->v; i++) {
	gchar *row[1];
        if (hidden_var(i, datainfo)) continue;
	row[0] = datainfo->varname[i];
        gtk_clist_append(GTK_CLIST(varbox), row);
    }
    gtk_clist_set_column_width (GTK_CLIST(varbox), 0, 80);
    gtk_widget_show(varbox);       
    gtk_box_pack_start(GTK_BOX(hbox), varbox, FALSE, TRUE, 0);

    vbox = gtk_vbox_new(FALSE, 1);

    /* this vbox will contain hboxes for dep var and indep vars */

    box = gtk_hbox_new(FALSE, 1);

    indepvarbox = gtk_clist_new(1);
    gtk_clist_clear(GTK_CLIST(indepvarbox));
    for (i=0; i<5; i++) {
	gchar *row[1];
	row[0] = " ";
        gtk_clist_append(GTK_CLIST(indepvarbox), row);
    }
    gtk_clist_set_column_width (GTK_CLIST(indepvarbox), 0, 80);
    gtk_widget_show(indepvarbox); 
    gtk_box_pack_start(GTK_BOX(vbox), indepvarbox, FALSE, TRUE, 0);
 
    gtk_widget_show(vbox);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, TRUE, 0);    

    gtk_widget_show(hbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, TRUE, 0);

    /* buttons: Estimate, Cancel, Help */
    tempwid = gtk_button_new_with_label (_("Estimate"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(dummy_call), NULL);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label(_("Cancel"));
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
                       GTK_SIGNAL_FUNC(delete_widget), dlg);
    gtk_widget_show(tempwid);

    tempwid = gtk_button_new_with_label(_("Help"));
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(context_help), 
			GINT_TO_POINTER(0));
    gtk_widget_show (tempwid);

    gtk_widget_show (dlg);
}
