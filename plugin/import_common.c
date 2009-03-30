/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

static void invalid_varname (PRN *prn)
{
    pputs(prn, gretl_errmsg_get());
    pputs(prn, _("\nPlease rename this variable and try again"));
}

#ifndef ODS_IMPORTER

#ifdef EXCEL_IMPORTER
# define cell_string(i,j) ((rows[i].cells != NULL)? rows[i].cells[j] : NULL)
#else
# define cell_string(i,j) (labels[i])
#endif

static int 
importer_dates_check (int row_offset, int col_offset, 
		      BookFlag *pflags, char **labels, 
		      DATAINFO *newinfo, PRN *prn, int *err)
{
    int d, t;
    char dstr[12];
    char *s;
    int ret = 0;

    for (t=0; t<newinfo->n; t++) {
	s = cell_string(t + row_offset, col_offset);
	if (s == NULL || *s == '\0') {
	    fprintf(stderr, "importer_dates_check: got blank label\n");
	    return 0;
	}
    }

    *err = dataset_allocate_obs_markers(newinfo);
    if (*err) {
	return 0;
    }

    for (t=0; t<newinfo->n && !*err; t++) {
	s = cell_string(t + row_offset, col_offset);
	if (*s == '"' || *s == '\'') s++;
	if (*pflags & BOOK_NUMERIC_DATES) {
	    if (sscanf(s, "%d", &d)) {
		MS_excel_date_string(dstr, d, 0, *pflags & BOOK_DATE_BASE_1904);
		s = dstr;
	    } else {
		pprintf(prn, "Bad date on row %d: '%s'\n", t+1, s);
		*err = E_DATA;
	    }
	}
	strncat(newinfo->S[t], s, OBSLEN - 1);
    }

    if (!*err) {
	int reversed = 0;

	ret = test_markers_for_dates(NULL, newinfo, &reversed, NULL, prn);
	if (reversed) {
	    *pflags |= BOOK_DATA_REVERSED;
	}
    }

    if (newinfo->markers != DAILY_DATE_STRINGS) {
	dataset_destroy_obs_markers(newinfo);
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
    free(book->targname);
    free(book->byte_offsets);
    free(book->xf_list);
}

static int wbook_check_params (wbook *book)
{
    if (book->targname != NULL) {
	int i;

	book->selected = -1;
	for (i=0; i<book->nsheets; i++) {
	    if (!strcmp(book->targname, book->sheetnames[i])) {
		book->selected = i;
	    }
	}
	if (book->selected < 0) {
	    gretl_errmsg_sprintf("\"%s\": no such sheet", book->targname);
	    return E_DATA;
	}
    }

    if (book->selected < 0 || book->selected >= book->nsheets) {
	return E_DATA;
    } else if (book->col_offset < 0 || book->row_offset < 0) {
	return E_DATA;
    } else {
	return 0;
    }
}

static void wbook_record_params (wbook *book, int *list)
{
    if (list != NULL && list[0] == 3) {
	list[1] = book->selected + 1;
	list[2] = book->col_offset;
	list[3] = book->row_offset;
    }
}

#endif /* !ODS_IMPORTER */

static void wbook_init (wbook *book, const int *list, char *sheetname)
{
    book->version = 0;
    book->nsheets = 0;
    book->col_offset = book->row_offset = 0;
    book->targname = NULL;
    book->sheetnames = NULL;
    book->byte_offsets = NULL;
    book->selected = 0;
    book->flags = 0;
    book->xf_list = NULL;
    book->get_min_offset = NULL;
    book->data = NULL;

    if (sheetname != NULL && *sheetname != '\0') {
	book->targname = gretl_strdup(sheetname);
	tailstrip(book->targname);
    }

    if (list != NULL && list[0] == 3) {
	if (book->targname == NULL && list[1] > 0) {
	    book->selected = list[1] - 1;
	}
	book->col_offset = list[2];
	book->row_offset = list[3];
    }
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

#ifdef EXCEL_IMPORTER
# ifndef WIN32

void infobox (const char *template, ...)
{
    GtkWidget *dialog;
    char msg[MAXLEN];
    va_list args;

    va_start(args, template);
    vsprintf(msg, template, args);
    va_end(args);

    dialog = gtk_message_dialog_new (NULL, 
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     GTK_MESSAGE_INFO,
				     GTK_BUTTONS_CLOSE,
				     msg);
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
}

# endif /* !WIN32 */

static 
void debug_callback (GtkWidget *w, wbook *book)
{
    static int done;

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	book_set_debug(book);
    }

    if (book_debugging(book) && !done) {
# ifdef WIN32
	gchar *msg;

	make_debug_fname();
	msg = g_strdup_printf(I_("Sending debugging output to %s"),
			      debug_fname);
	MessageBox(NULL, msg, "gretl", MB_OK | MB_ICONINFORMATION);
	g_free(msg);	
# else
	infobox(_("Sending debugging output to %s"), "stderr");
# endif
	done = 1;
    }
}

#endif /* EXCEL_IMPORTER */

static int book_get_min_offset (wbook *book, int k)
{
    if (book->get_min_offset != NULL) {
	return book->get_min_offset(book, k);
    } else {
	return 1;
    }
}

static
void wsheet_menu_select_row (GtkTreeSelection *selection, wbook *book)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    GtkTreePath *path;
    gint *idx;

    gtk_tree_selection_get_selected(selection, &model, &iter);
    path = gtk_tree_model_get_path(model, &iter);
    idx = gtk_tree_path_get_indices(path);

    if (book->selected != idx[0]) {
	int offmin, offcurr;

	book->selected = idx[0];

	offmin = book_get_min_offset(book, COL_OFFSET);
	offcurr = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(book->colspin));
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(book->colspin), offmin, 256);
	if (offcurr != offmin) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(book->colspin),
				      offmin);
	}

	offmin = book_get_min_offset(book, ROW_OFFSET);
	offcurr = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(book->rowspin));
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(book->rowspin), offmin, 256);
	if (offcurr != offmin) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(book->rowspin),
				      offmin);
	}	
    }
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

    gtk_tree_model_get_iter_first(model, &iter);
    gtk_tree_selection_select_iter(gtk_tree_view_get_selection(view), 
				   &iter);
}

static 
void wsheet_menu_cancel (GtkWidget *w, wbook *book)
{
    book->selected = -1;
}

static 
void wbook_set_col_offset (GtkWidget *w, wbook *book)
{
    book->col_offset = gtk_spin_button_get_value_as_int
	(GTK_SPIN_BUTTON(book->colspin)) - 1;
}

static 
void wbook_set_row_offset (GtkWidget *w, wbook *book)
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
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
    g_signal_connect(G_OBJECT(select), "changed",
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
    int offmin;

    w = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(w), _("gretl: spreadsheet import"));

    g_signal_connect_after(G_OBJECT(w), "delete_event",
			   G_CALLBACK(wsheet_menu_cancel), book);
    g_signal_connect(G_OBJECT(w), "destroy",  
		     G_CALLBACK(gtk_main_quit), NULL);

    vbox = GTK_DIALOG(w)->vbox;
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);

    /* selection of starting column and row */
    label = gtk_label_new(_("Start import at:"));
    gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    /* starting column spinner */
    tmp = gtk_label_new(_("column:"));
    offmin = book->col_offset + 1;
    c_adj = gtk_adjustment_new(offmin, offmin, 256, 1, 1, 0);
    book->colspin = gtk_spin_button_new(GTK_ADJUSTMENT(c_adj), 1, 0);
    g_signal_connect(c_adj, "value_changed",
		     G_CALLBACK(wbook_set_col_offset), book);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->colspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), book->colspin, FALSE, FALSE, 5);

    /* starting row spinner */
    tmp = gtk_label_new(_("row:"));
    offmin = book->row_offset + 1;
    r_adj = gtk_adjustment_new(offmin, offmin, 256, 1, 1, 0);
    book->rowspin = gtk_spin_button_new(GTK_ADJUSTMENT(r_adj), 1, 0);
    g_signal_connect(r_adj, "value_changed",
		     G_CALLBACK(wbook_set_row_offset), book);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->rowspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), book->rowspin, FALSE, FALSE, 5);

    /* column label feedback */
    hbox = gtk_hbox_new (FALSE, 5);
    gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, FALSE, 0);
    label = gtk_label_new("(A)");
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 5);
    g_signal_connect(GTK_EDITABLE(book->colspin), "changed",
		     G_CALLBACK(colspin_changed), label);

    /* choose the worksheet (if applicable) */
    if (multisheet) {
	add_sheets_list(vbox, book);
    }

#ifdef EXCEL_IMPORTER
    tmp = gtk_check_button_new_with_label(_("Produce debugging output"));
    g_signal_connect(G_OBJECT(tmp), "toggled", G_CALLBACK(debug_callback), 
		     book);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);
#endif

    hbox = GTK_DIALOG(w)->action_area;
    gtk_button_box_set_layout(GTK_BUTTON_BOX(hbox), GTK_BUTTONBOX_END);
    gtk_button_box_set_spacing(GTK_BUTTON_BOX(hbox), 10);

    /* Cancel button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    g_signal_connect(G_OBJECT (tmp), "clicked", 
		     G_CALLBACK(wsheet_menu_cancel), book);
    g_signal_connect_swapped(G_OBJECT (tmp), "clicked", 
			     G_CALLBACK (gtk_widget_destroy), 
			     G_OBJECT (w));

    /* OK button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    g_signal_connect_swapped(G_OBJECT (tmp), "clicked", 
			     G_CALLBACK (gtk_widget_destroy), 
			     G_OBJECT (w));
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_widget_grab_default(tmp);

    gtk_entry_set_activates_default(GTK_ENTRY(book->colspin), TRUE);
    gtk_entry_set_activates_default(GTK_ENTRY(book->rowspin), TRUE);

    gtk_widget_show_all(w);

    gtk_window_set_modal(GTK_WINDOW(w), TRUE);

    gtk_main();
}
