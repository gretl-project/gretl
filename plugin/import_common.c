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

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# include "gtk_compat.h"
#endif

#if GTK_MAJOR_VERSION >= 3
# include <gdk/gdkkeysyms-compat.h>
#else
# include <gdk/gdkkeysyms.h>
#endif

#include "gretl_zip.h"

static int check_imported_varname (char *vname, int row, int col,
				   PRN *prn)
{
    int err;

    gretl_charsub(vname, ' ', '_');

    err = check_varname(vname);

    if (err == VARNAME_BADCHAR) {
	/* the first byte was OK: try overwriting illegal 
	   characters with '_' 
	*/
	char *s = vname + 1;

	while (*s) {
	    if (!(isalpha((unsigned char) *s))  
		&& !(isdigit((unsigned char) *s))
		&& *s != '_') {
		*s = '_';
	    }
	    s++;
	}

	err = 0;
    }

    if (err) {
	if (row >= 0 && col >= 0) {
	    pprintf(prn, _("At row %d, column %d:\n"), row, col);
	}
	pputs(prn, gretl_errmsg_get());
	pputs(prn, _("\nPlease rename this variable and try again"));
    }

    return (err)? E_DATA : 0;
}

#ifndef EXCEL_IMPORTER /* FIXME? */

static void import_ts_check (DATASET *dset)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    int reversed = 0;
    int mpd = -1;

    mpd = test_markers_for_dates(dset, &reversed, NULL, prn);

    if (mpd > 0) {
	pputs(prn, _("taking date information from row labels\n\n"));
	if (dset->markers != DAILY_DATE_STRINGS) {
	    dataset_destroy_obs_markers(dset);
	}
	if (reversed) {
	    reverse_data(dset, prn);
	}
    } 

#if ODEBUG
    fprintf(stderr, "xinfo->dset->pd = %d\n", xinfo->dset->pd);
#endif

    if (dset->pd != 1 || strcmp(dset->stobs, "1")) { 
        dset->structure = TIME_SERIES;
    }

    gretl_print_destroy(prn);
}

#endif /* !EXCEL_IMPORTER */

#if defined(ODS_IMPORTER) || defined(XLSX_IMPORTER)

static char *get_absolute_path (const char *fname)
{
    char buf[FILENAME_MAX];
    char *s, *ret = NULL;

    s = getcwd(buf, sizeof buf - strlen(fname) - 2);
    if (s != NULL) {
	ret = g_strdup_printf("%s/%s", s, fname);
    }

    return ret;
}

static void remove_temp_dir (char *dname)
{
    const char *udir = gretl_dotdir();

# ifdef G_OS_WIN32
    char *fullpath = g_strdup_printf("%s%s", udir, dname);

    gretl_chdir(udir);
    win32_delete_dir(fullpath);
    g_free(fullpath);
# else
    gretl_chdir(udir);
    gretl_deltree(dname);
# endif
}

# ifdef G_OS_WIN32

static int gretl_make_tempdir (char *dname)
{
    strcpy(dname, ".gretl-ssheet-XXXXXX");
    mktemp(dname);

    if (*dname == '\0') {
	return E_FOPEN;
    } else {
	return gretl_mkdir(dname);
    }
}

# else

static int gretl_make_tempdir (char *dname)
{
    char *s;
    int err = 0;

    strcpy(dname, ".gretl-ssheet-XXXXXX");
    s = mkdtemp(dname);

    if (s == NULL) {
	gretl_errmsg_set_from_errno("gretl_make_tempdir");
	err = E_FOPEN;
    } 

    return err;
}

# endif /* G_OS_WIN32 or not */

/* For ODS and XLSX: unzip the target file in the user's
   "dotdir". On successful completion @dname holds the
   name of the temporary subdirectory, in the dotdir, 
   holding the contents of the zipfile.
*/

static int open_import_zipfile (const char *fname, char *dname,
				PRN *prn)
{
    const char *udir = gretl_dotdir();
    const char *real_fname = fname;
    char *recoded_fname = NULL;
    char *abspath = NULL;
    GError *gerr = NULL;
    FILE *fp;
    int err = 0;

    errno = 0;
    *dname = '\0';

    /* In case @fname is in UTF-8 but we're on MS Windows (or
       conceivably on an old non-UTF-8 *nix system), grab the
       appropriately recoded filename to pass into the
       zipfile apparatus, since in that context a filename
       that works with plain system stdio is expected.
    */

    fp = gretl_fopen_with_recode(fname, "r", &recoded_fname);
    if (fp == NULL) {
	return E_FOPEN;
    }

    fclose(fp);

    if (recoded_fname != NULL) {
	real_fname = recoded_fname;
    }

    /* by doing chdir, we may lose track of the file if 
       its path is relative */
    if (!g_path_is_absolute(real_fname)) {
	abspath = get_absolute_path(real_fname);
	if (abspath != NULL) {
	    real_fname = abspath;
	}
    }

    /* cd to user dir and make temporary dir */
    if (gretl_chdir(udir)) {
	gretl_errmsg_set_from_errno(udir);
	err = E_FOPEN;
    } else {
	err = gretl_make_tempdir(dname);
	if (!err) {
	    err = gretl_chdir(dname);
	    if (err) {
		gretl_remove(dname);
	    }
	}
    }
    
    if (!err) {
	err = gretl_unzip_file(real_fname, &gerr);
	if (gerr != NULL) {
	    pprintf(prn, "gretl_unzip_file: '%s'\n", gerr->message);
	    g_error_free(gerr);
	} else if (err) {
	    pprintf(prn, "gretl_unzip_file: err = %d\n", err);
	}
    }

    free(abspath);
    free(recoded_fname);

    return err;
}

/* check for spurious empty columns at the right of the sheet */

static int import_prune_columns (DATASET *dset)
{
    int allmiss = 1, ndel = 0;
    int i, t, err = 0;

    for (i=dset->v-1; i>0 && allmiss; i--) {
	for (t=0; t<dset->n; t++) {
	    if (!na(dset->Z[i][t])) {
		allmiss = 0;
		break;
	    }
	}
	if (allmiss) ndel++;
    }

    if (ndel == dset->v - 1) {
	gretl_errmsg_set(_("No numeric data were found"));
	err = E_DATA;
    } else if (ndel > 0) {
	fprintf(stderr, "Sheet has %d trailing empty variables\n", ndel);
	err = dataset_drop_last_variables(dset, ndel);
    }

    return err;
}

#else /* !ODS, !XLSX */

static int worksheet_start_dataset (DATASET *newinfo)
{
    if (newinfo->v == 1) {
	/* only the constant is present! */
	gretl_errmsg_set(_("No numeric data were found"));
	return E_DATA;
    } else {
	/* create import dataset */
	return start_new_Z(newinfo, 0);
    }
}

static int 
importer_dates_check (char **labels, BookFlag *pflags,
		      DATASET *newset, PRN *prn, 
		      int *err)
{
    int d, t;
    char dstr[12];
    char *s;
    int ret = 0;

    for (t=0; t<newset->n; t++) {
	s = labels[t];
	if (s == NULL || *s == '\0') {
	    fprintf(stderr, "importer_dates_check: got blank label\n");
	    return 0;
	}
    }

    *err = dataset_allocate_obs_markers(newset);
    if (*err) {
	return 0;
    }

    for (t=0; t<newset->n && !*err; t++) {
	s = labels[t];
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
	strncat(newset->S[t], s, OBSLEN - 1);
    }

    if (!*err) {
	int reversed = 0;

	ret = test_markers_for_dates(newset, &reversed, NULL, prn);
	if (reversed) {
	    *pflags |= BOOK_DATA_REVERSED;
	}
    }

    if (newset->markers != DAILY_DATE_STRINGS) {
	dataset_destroy_obs_markers(newset);
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

/* @list may contain sheet number, row and/or column offset;
   @sheetname may contain the name of a specific sheet; but
   both may be NULL
*/

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

    dialog = gtk_message_dialog_new(NULL, 
				    GTK_DIALOG_DESTROY_WITH_PARENT,
				    GTK_MESSAGE_INFO,
				    GTK_BUTTONS_CLOSE,
				    "%s",
				    msg);
    gtk_dialog_run(GTK_DIALOG (dialog));
    gtk_widget_destroy(dialog);
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
    g_object_unref(G_OBJECT(store));

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 0, NULL);
    column = gtk_tree_view_column_new_with_attributes(NULL,
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
    gtk_box_pack_start(GTK_BOX(vbox), hsep, FALSE, FALSE, 5);

    label = gtk_label_new(_("Sheet to import:"));
    gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 5);

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, 5);

    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				   GTK_POLICY_AUTOMATIC,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(sw), view); 
}

static void make_wmenu_modal (GtkWidget *w, gpointer p)
{
    gtk_window_set_modal(GTK_WINDOW(w), TRUE);
}

gboolean esc_cancels (GtkWidget *w, GdkEventKey *key, wbook *book)
{
    if (key->keyval == GDK_Escape) {
	if (book != NULL) {
	    book->selected = -1;
	}
        gtk_widget_destroy(w);
	return TRUE;
    } else {
	return FALSE;
    }
}

static void wsheet_menu (wbook *book, int multisheet)
{
    GtkWidget *w, *tmp, *label;
    GtkWidget *vbox, *hbox;
    GtkAdjustment *c_adj, *r_adj;
    int offmin;

    w = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(w), _("gretl: spreadsheet import"));

    g_signal_connect_after(G_OBJECT(w), "delete_event",
			   G_CALLBACK(wsheet_menu_cancel), book);
    g_signal_connect(G_OBJECT(w), "destroy",  
		     G_CALLBACK(gtk_main_quit), NULL);
    g_signal_connect(G_OBJECT(w), "realize",
		     G_CALLBACK(make_wmenu_modal), NULL);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(w));
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);

    /* selection of starting column and row */
    label = gtk_label_new(_("Start import at:"));
    gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* starting column spinner */
    tmp = gtk_label_new(_("column:"));
    offmin = book->col_offset + 1;
    c_adj = (GtkAdjustment *) gtk_adjustment_new(offmin, offmin, 256, 1, 1, 0);
    book->colspin = gtk_spin_button_new(c_adj, 1, 0);
    g_signal_connect(c_adj, "value_changed",
		     G_CALLBACK(wbook_set_col_offset), book);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->colspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), book->colspin, FALSE, FALSE, 5);

    /* starting row spinner */
    tmp = gtk_label_new(_("row:"));
    offmin = book->row_offset + 1;
    r_adj = (GtkAdjustment *) gtk_adjustment_new(offmin, offmin, 256, 1, 1, 0);
    book->rowspin = gtk_spin_button_new(r_adj, 1, 0);
    g_signal_connect(r_adj, "value_changed",
		     G_CALLBACK(wbook_set_row_offset), book);
    gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(book->rowspin),
				      GTK_UPDATE_IF_VALID);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), book->rowspin, FALSE, FALSE, 5);

    /* column label feedback */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    label = gtk_label_new("(A)");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
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
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
#endif

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(w));
    gtk_button_box_set_layout(GTK_BUTTON_BOX(hbox), GTK_BUTTONBOX_END);
    gtk_box_set_spacing(GTK_BOX(hbox), 10);

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
    gtk_widget_set_can_default(tmp, TRUE);
    gtk_widget_grab_default(tmp);

    g_signal_connect(G_OBJECT(w), "key-press-event", 
		     G_CALLBACK(esc_cancels), book);

    gtk_entry_set_activates_default(GTK_ENTRY(book->colspin), TRUE);
    gtk_entry_set_activates_default(GTK_ENTRY(book->rowspin), TRUE);

    gtk_widget_show_all(w);
    gtk_main();
}

