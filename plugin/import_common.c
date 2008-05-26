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

static int label_is_date (char *str, int *submax)
{
    int len = strlen(str);
    int i, d, sub = 0, pd = 0;
    char *p;
    double dd;

#if 0
    fprintf(stderr, "label_is_date: looking at '%s'\n", str);
#endif

    for (i=0; i<len; i++) {
	if (str[i] == ':' || str[i] == 'Q') {
	    str[i] = '.';
	    break;
	}
    }

    p = strchr(str, '.');

    if (len == 4 && sscanf(str, "%d", &d) && d > 0 && d < 3000) {
	pd = 1;
    } else if (p != NULL) {
	if (len == 6) {
	    if (sscanf(str, "%lf", &dd) && dd > 0 && dd < 3000) { 
		sub = atoi(p + 1);
		if (sub > 0 && sub < 5) {
		    pd = 4;
		}
	    }
	} else if (len == 7) {
	    if (sscanf(str, "%lf", &dd) && dd > 0 && dd < 3000) {
		sub = atoi(p + 1);
		if (sub > 0 && sub < 13) {
		    pd = 12;
		}
	    }
	}
	if (sub > *submax) {
	    *submax = sub;
	}
    }

    return pd;
}

static int pd_from_dmult (double dm)
{
    int pd = 0;

    if (dm > 250.0) {
	pd = 1;
    } else if (dm > 90.0) {
	pd = 4;
    } else if (dm > 29.0) {
	pd = 12;
    } else if (dm > 6.0) {
	pd = 52;
    } else if (dm > 1.4) {
	pd = 5;
    } else if (dm > 1.16) {
	pd = 6;
    } else if (dm > 0.9) {
	pd = 7;
    }

    return pd;
}

#define dmax(f) ((f == 1)? 366 : (f == 4)? 92 : 31)

static int 
calendar_missing_obs (double d1, double d0, int pd, BookFlag flags)
{
    int diff = d1 - d0;
    int mc = 0;

    if (pd == 52) {
	if (diff > 7) {
	    mc = (diff / 7) - 1;
	}
    } else if (pd == 7) {
	if (diff > 1) {
	    mc = diff - 1;
	}
    } else if (pd == 1 || pd == 4 || pd == 12) {
	if (diff > dmax(pd)) {
	    double xmc = diff / (365.0 / pd);

	    mc = floor(xmc - .5);
	}
    } else if ((pd == 5 || pd == 6) && diff > 1) {
	char dstr[12];
	int d, wday;

	mc = diff - 1;

	for (d=d0+1; d<d1; d++) {
	    MS_excel_date_string(dstr, d, 0, flags & BOOK_DATE_BASE_1904);
	    wday = get_day_of_week(dstr);
	    if (wday == 0) {
		/* Sunday: should be missing */
		mc--;
	    } else if (pd == 5 && wday == 6) {
		/* Saturday: missing if 5-day data */
		mc--;
	    }
	}
    }

    return mc;
}

#ifdef EXCEL_IMPORTER
# define cell_val(t,o) (rows[t].cells[o])
#else
# define cell_val(t,o) (labels[t])
#endif

static int pd_from_numeric_dates (int nrows, int row_offset, int col_offset,
				  char **labels, wbook *book)
{
    char dstr[12];
    int tstart = 1 + row_offset;
    int ndays, nobs = nrows - tstart;
    int mc, t, pd = 0;
    int d1, d0, dn;
    double dmult;
    char *test;

    fprintf(stderr, "check for consistent numeric dates in col %d (nobs = %d)\n", 
	    col_offset, nobs);

    /* look at first date */
    test = cell_val(tstart, col_offset);
    if (sscanf(test, "%d", &d0)) {
	MS_excel_date_string(dstr, d0, 0, book_base_1904(book));
	fprintf(stderr, "numeric date on row %d = %d (%s)\n", tstart, 
		d0, dstr);
    } else {
	fprintf(stderr, "failed to read starting\n");
	return 0;
    }

    /* look at last date */
    test = cell_val(nrows - 1, col_offset);
    if (sscanf(test, "%d", &dn)) {
	MS_excel_date_string(dstr, dn, 0, book_base_1904(book));
	fprintf(stderr, "numeric date on row %d = %d (%s)\n", nrows - 1, 
		dn, dstr);
    } else {
	fprintf(stderr, "failed to read ending date\n");
	return 0;
    }

    /* compare number of obs and calendar span of data */
    ndays = dn - d0 + 1;
    dmult = (double) ndays / nobs;
    fprintf(stderr, "Calendar interval = %d days\n", ndays);
    fprintf(stderr, "Calendar days per observation = %g\n", dmult);
    pd = pd_from_dmult(dmult);

    if (pd > 0) {
	fprintf(stderr, "provisional data frequency = %d\n", pd);
    } else {
	fputs("Can't make sense of this\n", stderr);
	return 0;
    }

    book->totmiss = 0;

    for (t=tstart; t<nrows; t++) {
	test = cell_val(t, col_offset);

	if (sscanf(test, "%d", &d1) != 1) {
	    fprintf(stderr, "Problem: blank cell at row %d\n", t + 1);
	    return 0;
	}

	if (t > tstart) {
	    mc = calendar_missing_obs(d1, d0, pd, book->flags);
	    if (mc > 0) {
		fprintf(stderr, "row %d: calendar gap = %d, %d values missing?\n", 
			t, d1 - d0, mc);
		book->totmiss += mc;
	    }
	}

	d0 = d1;
    }

    if (book->totmiss > 0) {
	int i, s = 0;

	fprintf(stderr, "Total missing values = %d\n", book->totmiss);
	book->missmask = calloc(nobs + book->totmiss, 1);

	if (book->missmask == NULL) {
	    fprintf(stderr, "Out of memory allocating missing obs mask\n");
	    return 0;
	}

	for (t=tstart; t<nrows; t++) {
	    test = cell_val(t, col_offset);
	    sscanf(test, "%d", &d1);
	    if (t > tstart) {
		mc = calendar_missing_obs(d1, d0, pd, book->flags);
		for (i=0; i<mc; i++) {
		    book->missmask[s++] = 1;
		}
	    }
	    d0 = d1;
	    s++;
	}
    }

    fprintf(stderr, "Setting data frequency = %d\n", pd);

    return pd;
}

#if 0 /* not yet */

static int 
new_consistent_date_labels (int nrows, int row_offset, int col_offset, 
			    char **labels, DATAINFO *newinfo, 
			    PRN *prn, int *err)
{
    int i, t, tstart = 1 + row_offset;
    char *s;
    int ret = 0;

    for (t=tstart; t<nrows; t++) {
	s = cell_val(t, col_offset);
	if (*s == '\0') {
	    return 0;
	}
    }

    *err = dataset_allocate_obs_markers(newinfo);
    if (*err) {
	return 0;
    }

    i = 0;
    for (t=tstart; t<nrows; t++) {
	s = cell_val(t, col_offset);
	if (*s == '"' || *s == '\'') s++;
	strncat(newinfo->S[i++], s, OBSLEN - 1);
    }

    ret = test_markers_for_dates(NULL, newinfo, NULL, prn);

    dataset_destroy_obs_markers(newinfo);

    return ret;
}

#endif

static int 
consistent_date_labels (int nrows, int row_offset, int col_offset, 
			char **labels)
{
    int t, tstart = 1 + row_offset;
    int pd = 0, pdbak = 0;
    int submax = 0;
    double x, xbak = 0.0;
    char *test;

    fprintf(stderr, "testing for consistent date labels in col %d\n", 
	    col_offset);

    for (t=tstart; t<nrows; t++) {
	test = cell_val(t, col_offset);

	if (*test == '\0') {
	    fprintf(stderr, " no: blank cell at row %d\n", t + 1);
	    return 0;
	} else if (*test == '"' || *test == '\'') {
	    test++;
	}

	pd = label_is_date(test, &submax);
	if (pd == 0) {
	    fprintf(stderr, " no: label '%s' on row %d is not a valid date\n", 
		    test, t + 1);
	    return 0;
	}

	if (pd == 12 && t - tstart > 3 && submax == 4) {
	    /* appeared to be monthly but really quarterly? */
	    pd = pdbak = 4;
	}

	x = atof(test);

	if (t > tstart) {
	    if (pd != pdbak) {
		fprintf(stderr, " no: got inconsistent data frequencies %d and %d\n",
			pdbak, pd);
		return 0;
	    }
	    if (x <= xbak) {
		fprintf(stderr, " no: got %g <= %g\n", x, xbak);
		return 0;
	    }
	    if (pd == 1 && x - xbak != 1) {
		fprintf(stderr, " no: got %g - %g = %g\n", x, xbak, x - xbak);
		return 0;
	    }		
	}

	pdbak = pd;
	xbak = x;
    }

    fprintf(stderr, " yes: data frequency = %d\n", pd);

    return pd;
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
    free(book->missmask);
}

static int wbook_check_params (wbook *book)
{
    if (book->targname != NULL) {
	int i;

	for (i=0; i<book->nsheets; i++) {
	    if (!strcmp(book->targname, book->sheetnames[i])) {
		book->selected = i;
	    }
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
    book->totmiss = 0;
    book->missmask = NULL;
    book->get_min_offset = NULL;
    book->data = NULL;

    if (sheetname != NULL && *sheetname != '\0') {
	book->targname = gretl_strdup(sheetname);
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
    c_adj = gtk_adjustment_new(offmin, offmin, 256, 1, 1, 1);
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
    r_adj = gtk_adjustment_new(offmin, offmin, 256, 1, 1, 1);
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


