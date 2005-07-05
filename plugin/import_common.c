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

#define BUILDING_PLUGIN
#include "../gui2/dialogs.h"

#ifdef EXCEL_IMPORTER

static void set_all_missing (double **Z, DATAINFO *pdinfo)
{
    int i, t;

    for (i=1; i<pdinfo->v; i++) {
	for (t=0; t<pdinfo->n; t++) {
	    Z[i][t] = NADBL;
	}
    }
}

#endif

static void time_series_setup (const char *s, DATAINFO *newinfo, int pd,
			       int *text_cols, int *time_series, 
			       int *label_strings)
{
    if (*s == '"' || *s == '\'') s++;

    newinfo->pd = pd;
    newinfo->structure = TIME_SERIES;

    strcpy(newinfo->stobs, s);
    colonize_obs(newinfo->stobs);
    fprintf(stderr, "stobs='%s'\n", newinfo->stobs);

    newinfo->sd0 = get_date_x(newinfo->pd, newinfo->stobs);

    if (text_cols != NULL) {
	*text_cols = 1;
    }

    *time_series = 1;
    *label_strings = 0;
}

static void invalid_varname (PRN *prn)
{
    pputs(prn, get_gretl_errmsg());
    pputs(prn, _("\nPlease rename this variable and try again"));
}

static int label_is_date (char *str)
{
    int len = strlen(str);
    int i, d, pd = 0;
    double dd, sub;

    if (*str == '"' || *str == '\'') {
	str++;
	len--;
    }

#if 0
    fprintf(stderr, "label_is_date: looking at '%s'\n", str);
#endif

    for (i=0; i<len; i++) {
	if (str[i] == ':' || str[i] == 'Q') {
	    str[i] = '.';
	    break;
	}
    }

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

static int obs_column_heading (const char *label)
{
    int ret = 0;

    if (label == NULL) {
	ret = 1;
    } else {
#if 0
	fprintf(stderr, "obs_column_heading: looking at '%s'\n", label);
#endif
	if (*label == '"') {
	    label++;
	}
	if (*label == '\0') {
	    ret = 1;    
	} else {
	    gchar *test = g_strdup(label);

	    lower(test);
	    if (strncmp(test, "obs", 3) == 0 ||
		strcmp(test, "date") == 0 ||
		strcmp(test, "year") == 0) {
		ret = 1;
	    }
	    g_free(test);
	}
    }

    return ret;
}

static void wbook_print_info (wbook *book) 
{
    int i;

    fprintf(stderr, "Found %d sheet%s\n", book->nsheets,
	    (book->nsheets > 1)? "s" : "");
    
    for (i=0; i<book->nsheets; i++) {
	if (book->byte_offsets != NULL) {
	    fprintf(stderr, "%d: '%s' at offset %u\n", i, 
		    book->sheetnames[i], book->byte_offsets[i]);
	} else {
	    fprintf(stderr, "%d: '%s'\n", i, book->sheetnames[i]);
	}
    }
}

static void wbook_free (wbook *book)
{
    int i;

    for (i=0; i<book->nsheets; i++) {
	free(book->sheetnames[i]);
    }
    free(book->sheetnames);
    free(book->byte_offsets);
}

static void wbook_init (wbook *book, int source)
{
    book->source = source;
    book->version = 0;
    book->nsheets = 0;
    book->col_offset = book->row_offset = 0;
    book->sheetnames = NULL;
    book->byte_offsets = NULL;
    book->selected = 0;
    book->debug = 0;
}

static const char *column_label (int col)
{
    const char *abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    static char label[5];
    char a1, a2;

    if (col < 26) {
	a2 = abc[col];
	sprintf(label, "(%c)", a2);
    } else {
	a1 = abc[col / 26 - 1];
	a2 = abc[col % 26];
	sprintf(label, "(%c%c)", a1, a2);
    }

    return label;
}

static void colspin_changed (GtkEditable *ed, GtkWidget *w)
{
    const gchar *text = gtk_entry_get_text(GTK_ENTRY(ed));

    if (text != NULL && isdigit((unsigned char) *text)) {
	int col = atoi(text);

	if (col > 0 && col < 257) {
	    gtk_label_set_text(GTK_LABEL(w), column_label(col - 1));
	}
    }
}

static 
void debug_callback (GtkWidget *w, wbook *book)
{
    static int done;

    book->debug = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
    if (book->debug && !done) {
	gchar *msg;

#ifdef WIN32
	make_debug_fname();
	msg = g_strdup_printf(_("Sending debugging output to %s"),
			      debug_fname);
#else
	msg = g_strdup_printf(_("Sending debugging output to %s"),
			      "stderr");
#endif
	infobox(msg);
	g_free(msg);
	done = 1;
    }
}

#if GTK_MAJOR_VERSION >= 2

static
void wsheet_menu_select_row (GtkTreeSelection *selection, wbook *book)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    GtkTreePath *path;
    gint *idx;

    gtk_tree_selection_get_selected (selection, &model, &iter);
    path = gtk_tree_model_get_path (model, &iter);
    idx = gtk_tree_path_get_indices(path);
    book->selected = idx[0];
}

static 
void wsheet_menu_make_list (GtkTreeView *view, wbook *book)
{
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeIter iter;
    int i;

    gtk_list_store_clear(GTK_LIST_STORE(model));
    gtk_tree_model_get_iter_first(model, &iter);
    
    for (i=0; i<book->nsheets; i++) {
        gtk_list_store_append(GTK_LIST_STORE(model), &iter);
        gtk_list_store_set(GTK_LIST_STORE(model), &iter, 
			   0, book->sheetnames[i], -1);
    }

    gtk_tree_model_get_iter_first (model, &iter);
    gtk_tree_selection_select_iter (gtk_tree_view_get_selection (view), 
				    &iter);
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

static 
void add_sheets_list (GtkWidget *vbox, wbook *book)
{
    GtkWidget *label, *view, *sw, *hsep;
    GtkListStore *store;
    GtkTreeSelection *select;
    GtkCellRenderer *renderer; 
    GtkTreeViewColumn *column;

    store = gtk_list_store_new(1, G_TYPE_STRING);
    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref (G_OBJECT(store));

    renderer = gtk_cell_renderer_text_new ();
    g_object_set (renderer, "ypad", 0, NULL);
    column = gtk_tree_view_column_new_with_attributes (NULL,
                                                       renderer,
                                                       "text", 
                                                       0, NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);   
    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode (select, GTK_SELECTION_SINGLE);
    g_signal_connect (G_OBJECT(select), "changed",
		      G_CALLBACK(wsheet_menu_select_row),
		      book);

    wsheet_menu_make_list(GTK_TREE_VIEW(view), book);

    /* now set up the widgets */

    hsep = gtk_hseparator_new();
    gtk_container_add(GTK_CONTAINER(vbox), hsep);

    label = gtk_label_new(_("Sheet to import:"));
    gtk_container_add(GTK_CONTAINER(vbox), label);

    sw = gtk_scrolled_window_new (NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, 5);

    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (sw),
                                    GTK_POLICY_AUTOMATIC,
                                    GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW (sw),
                                         GTK_SHADOW_IN);
    gtk_container_add (GTK_CONTAINER(sw), view); 
}

static void wsheet_menu (wbook *book, int multisheet)
{
    GtkWidget *w, *tmp, *label;
    GtkWidget *vbox, *hbox;
    GtkObject *c_adj, *r_adj;

    w = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    gtk_window_set_title(GTK_WINDOW(w), _("gretl: spreadsheet import"));

    g_signal_connect_after(G_OBJECT(w), "delete_event",
			   G_CALLBACK(wsheet_menu_cancel), book);
    g_signal_connect(G_OBJECT(w), "destroy",  
		     G_CALLBACK(gtk_main_quit), NULL);

    vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 5);
    gtk_container_add(GTK_CONTAINER(w), vbox);

    /* selection of starting column and row */
    label = gtk_label_new(_("Start import at:"));
    gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, FALSE, 5);

    hbox = gtk_hbox_new (FALSE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, TRUE, TRUE, 5);

    /* starting column spinner */
    tmp = gtk_label_new(_("column:"));
    c_adj = gtk_adjustment_new(1, 1, 256, 1, 1, 1);
    book->colspin = gtk_spin_button_new (GTK_ADJUSTMENT(c_adj), 1, 0);
    g_signal_connect (c_adj, "value_changed",
		      G_CALLBACK (wbook_get_col_offset), book);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->colspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->colspin, FALSE, FALSE, 5);

    /* starting row spinner */
    tmp = gtk_label_new(_("row:"));
    r_adj = gtk_adjustment_new(1, 1, 256, 1, 1, 1);
    book->rowspin = gtk_spin_button_new (GTK_ADJUSTMENT(r_adj), 1, 0);
    g_signal_connect (r_adj, "value_changed",
		      G_CALLBACK (wbook_get_row_offset), book);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->rowspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->rowspin, FALSE, FALSE, 5);

    /* column label feedback */
    hbox = gtk_hbox_new (FALSE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, FALSE, 0);
    label = gtk_label_new("(A)");
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 5);
    g_signal_connect (GTK_EDITABLE(book->colspin), "changed",
		      G_CALLBACK (colspin_changed), label);

    /* choose the worksheet (if applicable) */
    if (multisheet) {
	add_sheets_list(vbox, book);
    }

    /* debugging option for XLS */
    if (book->source == WBOOK_XLS) {
	GtkWidget *chk;

	chk = gtk_check_button_new_with_label(_("Produce debugging output"));
	g_signal_connect(G_OBJECT(chk), "toggled", G_CALLBACK(debug_callback), 
			 book);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), FALSE);
	gtk_box_pack_start(GTK_BOX(vbox), chk, TRUE, TRUE, 5);
    }

    hbox = gtk_hbox_new (TRUE, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);

    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, 0);
    g_signal_connect_swapped (G_OBJECT (tmp), "clicked", 
			      G_CALLBACK (gtk_widget_destroy), 
			      G_OBJECT (w));
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_widget_grab_default(tmp);

    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, 0);
    g_signal_connect (G_OBJECT (tmp), "clicked", 
		      G_CALLBACK(wsheet_menu_cancel), book);
    g_signal_connect_swapped (G_OBJECT (tmp), "clicked", 
			      G_CALLBACK (gtk_widget_destroy), 
			      G_OBJECT (w));

    gtk_entry_set_activates_default(GTK_ENTRY(book->colspin), TRUE);
    gtk_entry_set_activates_default(GTK_ENTRY(book->rowspin), TRUE);

    gtk_widget_show_all(w);

    gtk_window_set_modal(GTK_WINDOW(w), TRUE);

    gtk_main();
}

#else /* GTK 1.2 */

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

static gint entry_activate (GtkWidget *w, GdkEventKey *key, gpointer p)
{
    GtkWidget *top = gtk_widget_get_toplevel(w);

    gtk_window_activate_default(GTK_WINDOW(top));

    return FALSE;
}

static void wsheet_menu (wbook *book, int multisheet)
{
    GtkWidget *w, *tmp, *frame;
    GtkWidget *vbox, *hbox, *list;
    GtkObject *c_adj, *r_adj;

    w = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    gtk_window_set_title(GTK_WINDOW(w), "gretl: spreadsheet import");

    gtk_signal_connect(GTK_OBJECT(w), "destroy",  
		       GTK_SIGNAL_FUNC(gtk_main_quit), NULL);

    vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 10);
    gtk_container_add(GTK_CONTAINER(w), vbox);    

    /* selection of starting column and row */
    tmp = gtk_label_new(_("Start import at:"));
    gtk_box_pack_start (GTK_BOX (vbox), tmp, FALSE, FALSE, 5);

    hbox = gtk_hbox_new (FALSE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, TRUE, TRUE, 5);

    /* starting column spinner */
    tmp = gtk_label_new("column:");
    c_adj = gtk_adjustment_new(1, 1, 256, 1, 1, 1);
    book->colspin = gtk_spin_button_new (GTK_ADJUSTMENT(c_adj), 1, 0);
    gtk_signal_connect (c_adj, "value_changed",
			GTK_SIGNAL_FUNC (wbook_get_col_offset), book);
    gtk_signal_connect(GTK_OBJECT(book->colspin), "activate", 
		       GTK_SIGNAL_FUNC(entry_activate), NULL);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->colspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->colspin, FALSE, FALSE, 5);

    /* starting row spinner */
    tmp = gtk_label_new("row:");
    r_adj = gtk_adjustment_new(1, 1, 256, 1, 1, 1);
    book->rowspin = gtk_spin_button_new (GTK_ADJUSTMENT(r_adj), 1, 0);
    gtk_signal_connect (r_adj, "value_changed",
			GTK_SIGNAL_FUNC (wbook_get_row_offset), book);
    gtk_signal_connect(GTK_OBJECT(book->rowspin), "activate", 
		       GTK_SIGNAL_FUNC(entry_activate), NULL);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->rowspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->rowspin, FALSE, FALSE, 5);

    /* column label feedback */
    hbox = gtk_hbox_new (FALSE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, FALSE, 0);
    tmp = gtk_label_new("(A)");
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_signal_connect (GTK_OBJECT(GTK_EDITABLE(book->colspin)), "changed",
			GTK_SIGNAL_FUNC (colspin_changed), tmp);

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

    /* debugging option for XLS */
    if (book->source == WBOOK_XLS) {
	GtkWidget *chk;

	chk = gtk_check_button_new_with_label(_("Produce debugging output"));
	gtk_signal_connect(GTK_OBJECT(chk), "toggled", 
			   GTK_SIGNAL_FUNC(debug_callback), 
			   book);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), FALSE);
	gtk_box_pack_start(GTK_BOX(vbox), chk, TRUE, TRUE, 5);
    }

    hbox = gtk_hbox_new (TRUE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, FALSE, 5);

    tmp = gtk_button_new_with_label("OK");
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));
    gtk_widget_grab_default (tmp);

    tmp = gtk_button_new_with_label("Cancel");
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC (wsheet_menu_cancel), book);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));

    gtk_widget_show_all(w);

    gtk_main();
}

#endif /* gtk 2.0 vs 1.2 */

