/*
 *  Copyright (c) by Allin Cottrell 2002
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

static void set_all_missing (double **Z, DATAINFO *pdinfo)
{
    int i, t;

    for (i=1; i<pdinfo->v; i++)
	for (t=0; t<pdinfo->n; t++)
	    Z[i][t] = NADBL;
}

static int label_is_date (char *str)
{
    size_t len = strlen(str);
    int i, d, pd = 0;
    double dd, sub;

    for (i=0; i<len; i++)
	if (str[i] == ':') str[i] = '.';
     
    if (len == 4 && sscanf(str, "%4d", &d) == 1 &&
	d > 0 && d < 3000) {
	pd = 1;
    }
    else if (len == 6 && sscanf(str, "%lf", &dd) == 1 &&
	dd > 0 && dd < 3000) { 
	sub = 10.0 * (dd - (int) dd);
	if (sub >= .999 && sub <= 4.001) pd = 4;
    }
    else if (len == 7 && sscanf(str, "%lf", &dd) == 1 &&
	dd > 0 && dd < 3000) {
	sub = 100.0 * (dd - (int) dd);
	if (sub >= .9999 && sub <= 12.0001) pd = 12;
    }
    return pd;
}

static int obs_column (char *label)
{
    fprintf(stderr, "obs_column(): test='%s'\n", label);

    if (*label == '\0') return 1;    

    lower(label);
    if (strcmp(label, "obs") == 0 ||
	strcmp(label, "date") == 0 ||
	strcmp(label, "year") == 0)
	return 1;

    return 0;
}

static void wbook_init (wbook *book)
{
    book->nsheets = 0;
    book->col_offset = book->row_offset = 0;
    book->sheetnames = NULL;
    book->selected = 0;
}

static
void wsheet_menu_select_row (GtkCList *clist, gint row, gint column, 
			     GdkEventButton *event, wbook *book) 
{
    book->selected = row;
}

static 
void wsheet_menu_make_list (GtkWidget *widget, wbook *book)
{
    gchar *row[1];
    int i;

    gtk_clist_clear(GTK_CLIST(widget));
    for (i=0; i<book->nsheets; i++) {
	row[0] = book->sheetnames[i];
        gtk_clist_append(GTK_CLIST (widget), row);
    }

    gtk_clist_select_row(GTK_CLIST(widget), 0, 0);  
}

static 
void wsheet_menu_cancel (GtkWidget *w, wbook *book)
{
    book->selected = -1;
}

static 
void wbook_get_col_offset (GtkWidget *w, wbook *book)
{
    book->col_offset = gtk_spin_button_get_value_as_int
	(GTK_SPIN_BUTTON(book->colspin)) - 1;
}

static 
void wbook_get_row_offset (GtkWidget *w, wbook *book)
{
    book->row_offset = gtk_spin_button_get_value_as_int
	(GTK_SPIN_BUTTON(book->rowspin)) - 1;
}

static void wsheet_menu (wbook *book, int multisheet)
{
    GtkWidget *w, *tmp, *frame;
    GtkWidget *vbox, *hbox, *list;
    GtkObject *adj;

    w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(w), "gretl: spreadsheet import");
    gtk_signal_connect(GTK_OBJECT(w), "destroy",  
		       GTK_SIGNAL_FUNC(gtk_main_quit), NULL);

    vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 10);

    /* choose starting column and row */
    frame = gtk_frame_new("Start import at");
    gtk_box_pack_start (GTK_BOX (vbox), frame, TRUE, TRUE, 5);

    hbox = gtk_hbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (hbox), 10);
    gtk_container_add(GTK_CONTAINER (frame), hbox);

    tmp = gtk_label_new("column:");
    adj = gtk_adjustment_new(1, 1, 5, 1, 1, 1);
    book->colspin = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    gtk_signal_connect (adj, "value_changed",
			GTK_SIGNAL_FUNC (wbook_get_col_offset), book);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->colspin, FALSE, FALSE, 5);

    tmp = gtk_label_new("row:");
    adj = gtk_adjustment_new(1, 1, 5, 1, 1, 1);
    book->rowspin = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    gtk_signal_connect (adj, "value_changed",
			GTK_SIGNAL_FUNC (wbook_get_row_offset), book);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->rowspin, FALSE, FALSE, 5);

    /* choose the worksheet (if applicable) */
    if (multisheet) {
	frame = gtk_frame_new("Sheet to import");
	gtk_box_pack_start (GTK_BOX (vbox), frame, TRUE, TRUE, 5);
	list = gtk_clist_new(1);
	gtk_signal_connect (GTK_OBJECT (list), "select_row", 
			    GTK_SIGNAL_FUNC (wsheet_menu_select_row), 
			    book);
	wsheet_menu_make_list(list, book);
	gtk_container_set_border_width (GTK_CONTAINER (list), 10);
	gtk_container_add(GTK_CONTAINER(frame), list);
    }

    hbox = gtk_hbox_new (TRUE, 5);
    tmp = gtk_button_new_with_label("OK");
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));
    gtk_widget_grab_default (tmp);

    tmp = gtk_button_new_with_label("Cancel");
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC (wsheet_menu_cancel), book);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));

    gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, FALSE, 5);

    gtk_container_add(GTK_CONTAINER(w), vbox);

    gtk_widget_show_all(w);
    gtk_main();
}

