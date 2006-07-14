/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

#include "gretl.h"

#include <errno.h>
#include <ctype.h>

#include <gtkextra/gtksheet.h>
#include <gtkextra/gtkitementry.h>

#include "ssheet.h"
#include "dlgutils.h"
#include "menustate.h"

#define DEFAULT_PRECISION 6
#define DEFAULT_SPACE 8

static GtkWidget *gretlsheet;
static GtkWidget *sheetwin;
static GtkWidget *locator;
static GtkWidget *topentry;
static GtkWidget *sheet_popup;

static char newvarname[VNAMELEN];
static char newobsmarker[OBSLEN];
static int numcolumns, numrows;
static int sheet_modified;
static int origvars;
static int orig_main_v;
static int sheetcode;

static void sheet_add_var_callback (void);
static void sheet_add_obs_callback (void);
static void sheet_insert_obs_callback (void);
static void sheet_clear (gpointer *data, guint all, GtkWidget *w);
static void get_data_from_sheet (void);
static void set_dataset_locked (gboolean s);

static GtkItemFactoryEntry sheet_items[] = {
    { N_("/_Observation"), NULL, NULL, 0, "<Branch>" },
    { N_("/Observation/_Append obs"), NULL, sheet_add_obs_callback, 0, NULL },
    { N_("/Observation/_Insert obs"), NULL, sheet_insert_obs_callback, 0, NULL },
    { N_("/_Variable"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/_Add"), NULL, sheet_add_var_callback, 0, NULL },
    { N_("/_Clear"), NULL, NULL, 0, "<Branch>" },
    { N_("/_Clear/_Selected cells"), NULL, sheet_clear, 0, NULL },
    { N_("/_Clear/_All data"), NULL, sheet_clear, 1, NULL }
};

enum {
    SHEET_AT_END,
    SHEET_AT_POINT
};

static void disable_obs_menu (GtkItemFactory *ifac)
{
    GtkWidget *w = gtk_item_factory_get_item(ifac, "/Observation");

    gtk_widget_set_sensitive(w, FALSE);
}

static int spreadsheet_hide (int i, const DATAINFO *pdinfo)
{
    int ret = 0;

    if (var_is_scalar(pdinfo, i)) {
	ret = 1;
    } else if (var_is_hidden(pdinfo, i)) {
	ret = 1;
    }

    return ret;
}

static int check_data_in_sheet (void)
{
    gint err, i, t, n = datainfo->n;
    gchar *celltext;
    char numstr[32];

    for (i=0; i<numcolumns; i++) { 
	for (t=0; t<n; t++) {
	    celltext = gtk_sheet_cell_get_text(GTK_SHEET(gretlsheet), t, i);
	    if (celltext != NULL) {
		*numstr = 0;
		strncat(numstr, celltext, 31);
		err = check_atof(numstr);
		if (err) {
		    errbox(get_gretl_errmsg());
		    return 1;
		}
	    }
	}
    }

    return 0;
}

static void real_add_new_var (GtkWidget *widget, dialog_t *ddata) 
{
    GtkSheet *sheet = GTK_SHEET(gretlsheet);     
    const char *buf = edit_dialog_get_text(ddata);

    if (buf == NULL || validate_varname(buf)) {
	return;
    }

    *newvarname = '\0';
    strncat(newvarname, buf, VNAMELEN - 1);

    close_dialog(ddata);

    if (*newvarname != '\0') {
	gtk_sheet_add_column(sheet, 1);
	gtk_sheet_column_button_add_label(sheet, numcolumns, newvarname);
	gtk_sheet_set_column_title(sheet, numcolumns, newvarname);
	gtk_sheet_set_active_cell(sheet, 0, numcolumns);
	gtk_sheet_moveto(sheet, 0, numcolumns, 0.0, 0.8);
	numcolumns++;
	sheet_modified = 1;
    }
}

static void real_add_named_obs (dialog_t *ddata, int flag) 
{
    GtkSheet *sheet = GTK_SHEET(gretlsheet);
    const char *buf = edit_dialog_get_text(ddata);

    if (buf == NULL) return;

    *newobsmarker = '\0';
    strncat(newobsmarker, buf, OBSLEN - 1);

    close_dialog(ddata);

    if (*newobsmarker != '\0') {
	int addrow;

	if (flag == SHEET_AT_POINT) {
	    addrow = sheet->active_cell.row;
	    gtk_sheet_insert_rows(sheet, addrow, 1);
	} else {
	    addrow = numrows;
	    gtk_sheet_add_row(sheet, 1);
	}

	gtk_sheet_row_button_add_label(sheet, addrow, newobsmarker);
	gtk_sheet_set_row_title(sheet, addrow, newobsmarker);
	gtk_sheet_set_active_cell(sheet, addrow, 0);
	gtk_sheet_moveto(sheet, addrow, 0, 0.8, 0.0);
	numrows++;
	sheet_modified = 1;
    }
}

static void insert_named_obs (GtkWidget *widget, dialog_t *ddata)
{
    real_add_named_obs(ddata, SHEET_AT_POINT);
}

static void append_named_obs (GtkWidget *widget, dialog_t *ddata)
{
    real_add_named_obs(ddata, SHEET_AT_END);
}

static void name_var_dialog (void) 
{
    edit_dialog (_("gretl: name variable"), 
		 _("Enter name for new variable\n"
		 "(max. 15 characters)"),
		 NULL, real_add_new_var, mdata, 
		 0, VARCLICK_NONE, NULL);
}

static void name_obs_dialog (int flag) 
{
    edit_dialog (_("gretl: case marker"), 
		 _("Enter case marker for new obs\n"
		 "(max. 15 characters)"),
		 NULL, 
		 (flag == SHEET_AT_POINT)? insert_named_obs : append_named_obs, 
		 mdata, 
		 0, VARCLICK_NONE, NULL);
}

static void sheet_add_var_callback (void)
{
    *newvarname = '\0';
    name_var_dialog();
}

static void sheet_add_obs_callback (void)
{
    GtkSheet *sheet = GTK_SHEET(gretlsheet);
    char rowlabel[OBSLEN];
    int i, n = 1;

    if (datainfo->markers) {
	name_obs_dialog(SHEET_AT_END);
	return;
    } else {
	n = add_obs_dialog(NULL, 1);
    }

    if (n <= 0) return;

    for (i=0; i<n; i++) {
	gtk_sheet_add_row(sheet, 1);
	get_full_obs_string(rowlabel, numrows + i, datainfo);
	gtk_sheet_row_button_add_label(sheet, numrows + i, rowlabel);
	gtk_sheet_set_row_title(sheet, numrows + i, rowlabel);
    }

    gtk_sheet_set_active_cell(sheet, numrows, 0);
    gtk_sheet_moveto(sheet, numrows, 0, 0.8, 0.0);
    numrows += n;
    sheet_modified = 1;
}

static void sheet_insert_obs_callback (void)
{
    GtkSheet *sheet = GTK_SHEET(gretlsheet);
    char rowlabel[10];
    int i;

    if (datainfo->markers) {
	name_obs_dialog(SHEET_AT_POINT);
	return;
    }

    gtk_sheet_insert_rows(sheet, sheet->active_cell.row, 1);
    numrows++;

    for (i=sheet->active_cell.row; i<numrows; i++) {
	get_full_obs_string(rowlabel, i, datainfo);
	gtk_sheet_row_button_add_label(sheet, i, rowlabel);
	gtk_sheet_set_row_title(sheet, i, rowlabel);
    }

    sheet_modified = 1;
}

static void sheet_clear (gpointer *data, guint all, GtkWidget *w)
{
    GtkSheet *sheet = GTK_SHEET(gretlsheet);

    if (!all && sheet->state != GTK_SHEET_NORMAL) { 
	gtk_sheet_range_clear(sheet, &sheet->range);
    } else {
	gtk_sheet_range_clear(sheet, NULL);
    }	
}

static gint popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("Add Variable")) == 0) { 
	sheet_add_var_callback();
    } else if (strcmp(item, _("Add Observation")) == 0) {
	sheet_add_obs_callback();
    } else if (strcmp(item, _("Insert Observation")) == 0) {
	sheet_insert_obs_callback();
    } else if (strcmp(item, _("Clear Cells")) == 0) {
	sheet_clear(NULL, 0, NULL);
    } 

    gtk_widget_destroy(sheet_popup);

    return TRUE;
}

static GtkWidget *build_menu (GtkWidget *sheet)
{
    static char *items[]={
        N_("Add Variable"),
        N_("Add Observation"),
        N_("Insert Observation"),
        N_("Clear Cells")
    };
    GtkWidget *menu;
    GtkWidget *item;
    int i, n_items = sizeof items / sizeof items[0];

    menu = gtk_menu_new();

    for (i=0; i<n_items; i++) {
        item = gtk_menu_item_new_with_label(_(items[i]));
        gtk_signal_connect(GTK_OBJECT(item), "activate",
                           (GtkSignalFunc) popup_activated,
                           _(items[i]));
        GTK_WIDGET_SET_FLAGS (item, GTK_SENSITIVE | GTK_CAN_FOCUS);
        gtk_widget_show(item);
        gtk_menu_append(GTK_MENU(menu), item);
    }

    return menu;
}

static gint do_sheet_popup (GtkWidget *widget, GdkEventButton *event, 
			    gpointer data)
{
    GdkModifierType mods;
    GtkWidget *sheet;

    sheet = GTK_WIDGET(widget);

    gdk_window_get_pointer (sheet->window, NULL, NULL, &mods);
    if (mods&GDK_BUTTON3_MASK) { 
        if (sheet_popup) g_free(sheet_popup);
        sheet_popup = build_menu(sheet);
        gtk_menu_popup(GTK_MENU(sheet_popup), NULL, NULL, NULL, NULL,
                       event->button, event->time);
    }

    return TRUE;
}

static void parse_numbers (GtkWidget *widget, gpointer data)
{
    GtkSheet *sheet = GTK_SHEET(widget);
    gchar *entrytext;
    char label[32];

    entrytext = gtk_entry_get_text(GTK_ENTRY(sheet->sheet_entry));
    if (entrytext == NULL) return;

    if (!strcmp(entrytext, "na") || !strcmp(entrytext, "NA")) {
	*label = 0;
    } else if (check_atof(entrytext)) {
	errbox(get_gretl_errmsg());
	*label = 0;
    } else {
	sprintf(label, "%.*f", DEFAULT_PRECISION, atof(entrytext));
    }

    gtk_sheet_set_cell(sheet, sheet->active_cell.row,
                       sheet->active_cell.col, GTK_JUSTIFY_RIGHT, label); 

    sheet_modified = 1;
}

static void show_sheet_entry (GtkWidget *widget, gpointer data)
{
    gchar *text;
    GtkSheet *sheet;
    GtkEntry *sheet_entry;

    if (!GTK_WIDGET_HAS_FOCUS(widget)) return;

    sheet = GTK_SHEET(gretlsheet);
    sheet_entry = GTK_ENTRY(sheet->sheet_entry);

    if ((text = gtk_entry_get_text(GTK_ENTRY(topentry)))) {
        gtk_entry_set_text(sheet_entry, text);
    }
}

static void activate_sheet_entry (GtkWidget *widget, gpointer data)
{
    GtkSheet *sheet = GTK_SHEET(gretlsheet);
    GtkEntry *sheet_entry;
    gint row, col;
    gint just = GTK_JUSTIFY_RIGHT;
  
    row = sheet->active_cell.row; 
    col = sheet->active_cell.col;

    sheet_entry = GTK_ENTRY(gtk_sheet_get_entry(sheet));

    if (GTK_IS_ITEM_ENTRY(sheet_entry))
        just = GTK_ITEM_ENTRY(sheet_entry)->justification;

    gtk_sheet_set_cell(sheet, row, col,
                       just, gtk_entry_get_text(sheet_entry));

}

static void show_entry (GtkWidget *widget, gpointer data)
{
    gchar *text; 
    GtkSheet *sheet = GTK_SHEET(gretlsheet);

    if (!GTK_WIDGET_HAS_FOCUS(widget)) return;

    if ((text = gtk_entry_get_text(GTK_ENTRY(sheet->sheet_entry)))) {
	chopstr(text);
	gtk_entry_set_text(GTK_ENTRY(topentry), text);
    }
}

static gint activate_sheet_cell (GtkWidget *w, gint row, gint column, 
				 gpointer data) 
{
    GtkSheet *sheet = GTK_SHEET(w);
    GtkEntry *sheet_entry;
    char cell[100];
    char *text;

    sheet_entry = GTK_ENTRY(gtk_sheet_get_entry(sheet));

    if (sheet->column[column].name && sheet->row[row].name) {
        sprintf(cell,"  %s:%s  ", 
                sheet->column[column].name, 
                sheet->row[row].name);
    } else {
        sprintf(cell,"  %s:%d  ", sheet->column[column].name, row);
    }

    gtk_label_set(GTK_LABEL(locator), cell);
    gtk_entry_set_max_length(GTK_ENTRY(topentry),
                             GTK_ENTRY(sheet_entry)->text_max_length);

    if ((text = gtk_entry_get_text(GTK_ENTRY(gtk_sheet_get_entry(sheet))))) {
        gtk_entry_set_text(GTK_ENTRY(topentry), text);
    } else {
        gtk_entry_set_text(GTK_ENTRY(topentry), "");
    }

    return TRUE;
}

static int var_added_since_ssheet_opened (int i, int main_v)
{
    return (i >= orig_main_v && i < main_v);
}

static void get_data_from_sheet (void)
{
    int oldv = datainfo->v;
    int newvars = numcolumns - origvars;
    int newobs = numrows - datainfo->n;
    int missobs = 0;

    gchar *celltext;
    GtkSheet *sheet = GTK_SHEET(gretlsheet);

    int i, colnum, t;

    if (check_data_in_sheet()) {
	return;
    }

    if (!sheet_modified) {
	infobox(_("No changes were made"));
	return;
    }

    if (newobs > 0) {
	if (dataset_add_observations(newobs, &Z, datainfo, OPT_A) ||
	    dataset_destroy_hidden_variables(&Z, datainfo, 0)) {
	    errbox(_("Failed to allocate memory for new data"));
	    return;
	}
    }

    if (newvars > 0) {
	if (dataset_add_series(newvars, &Z, datainfo)) {
	    errbox(_("Failed to allocate memory for new data"));
	    return;
	}
	for (i=0; i<newvars; i++) { 
	    colnum = origvars + i; 
	    strcpy(datainfo->varname[oldv + i], sheet->column[colnum].button.label);
	}
    }

    colnum = 0;

    for (i=1; i<datainfo->v; i++) {
	if (spreadsheet_hide(i, datainfo)) {
	    continue;
	}
	if (var_added_since_ssheet_opened(i, oldv)) {
	    continue;
	}
	for (t=0; t<datainfo->n; t++) {
	    celltext = gtk_sheet_cell_get_text(sheet, t, colnum);
	    if (celltext != NULL && *celltext != 0) {
		Z[i][t] = atof(celltext);
	    } else {
		Z[i][t] = NADBL;
		missobs = 1;
	    }
	}
	colnum++;
    }

    if (datainfo->markers && datainfo->S != NULL) {
	for (t=0; t<datainfo->n; t++) {
	    strcpy(datainfo->S[t], sheet->row[t].button.label); 
	}
    }

    data_status |= (GUI_DATA|MODIFIED_DATA);
    register_data(NULL, NULL, 0);

    if (missobs) {
	infobox(_("Warning: there were missing observations"));
    } 

    sheet_modified = 0;
}

static void add_data_to_sheet (GtkWidget *w)
{
    gchar label[32], rowlabel[10];
    GtkSheet *sheet = GTK_SHEET(w);
    int i, j, t;

    numrows = datainfo->n;
    numcolumns = j = 0;

    for (i=1; i<datainfo->v; i++) {
	if (spreadsheet_hide(i, datainfo)) {
	    continue;
	}
	gtk_sheet_column_button_add_label(sheet, j, datainfo->varname[i]);
	gtk_sheet_set_column_title(sheet, j, datainfo->varname[i]);
	j++;
    }

    numcolumns = origvars = j;

    for (t=0; t<datainfo->n; t++) {
	get_full_obs_string(rowlabel, t, datainfo);
	gtk_sheet_row_button_add_label(sheet, t, rowlabel);
	gtk_sheet_set_row_title(sheet, t, rowlabel);
    }

    j = 0;
    for (i=1; i<datainfo->v; i++) {
	if (spreadsheet_hide(i, datainfo)) {
	    continue;
	}
	for (t=0; t<datainfo->n; t++) {
	    if (!na(Z[i][t])) {
		sprintf(label, "%.*f", DEFAULT_PRECISION, Z[i][t]);
	    } else {
		*label = 0;
	    }
	    gtk_sheet_set_cell(sheet, t, j, GTK_JUSTIFY_RIGHT, label); 
	}
	j++;
    }

    orig_main_v = datainfo->v;
    
    sheet_modified = 0;
}

static void add_skel_to_sheet (GtkWidget *w)
{
    GtkSheet *sheet = GTK_SHEET(w);
    gchar rowlabel[10];
    int i, j, t;

    numrows = datainfo->n;
    numcolumns = j = 0;

    for (i=1; i<datainfo->v; i++) {
	gtk_sheet_column_button_add_label(sheet, j, datainfo->varname[i]);
	gtk_sheet_set_column_title(sheet, j, datainfo->varname[i]);
	j++;
    }

    numcolumns = origvars = j;

    for (t=0; t<datainfo->n; t++) {
	get_full_obs_string(rowlabel, t, datainfo);
	gtk_sheet_row_button_add_label(sheet, t, rowlabel);
	gtk_sheet_set_row_title(sheet, t, rowlabel);
    }

    j = 0;
    for (i=1; i<datainfo->v; i++) {
	for (t=0; t<datainfo->n; t++) {
	    gtk_sheet_set_cell(sheet, t, j, GTK_JUSTIFY_RIGHT, ""); 
	}
	j++;
    }

    sheet_modified = 0;
}

static void free_spreadsheet (GtkWidget *widget, gpointer data) 
{
    gretlsheet = NULL;
    set_dataset_locked(FALSE);
}

static void empty_dataset_guard (void)
{
    int t, empty = 0;

    if (datainfo->v == 2) {
	empty = 1;
	for (t=0; t<datainfo->n; t++) {
	    if (!na(Z[1][t])) {
		empty = 0;
		break;
	    }
	}
    }

    if (empty) {
	data_status |= (GUI_DATA | MODIFIED_DATA);
	register_data(NULL, NULL, 0);
	infobox(_("Warning: there were missing observations"));
    }
}

static void maybe_exit_sheet (GtkWidget *w, GtkWidget *swin)
{
    int resp;

    if (sheet_modified) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    get_data_from_sheet();
	} else if (resp == GRETL_CANCEL || resp == -1) {
	    return;
	}
    } 

    if (sheetcode == SHEET_NEW_DATASET) {
	empty_dataset_guard();
    }
  
    gtk_widget_destroy(swin);
}

static gint sheet_delete_event (GtkWidget *w, GdkEvent *event,
				gpointer p)
{
    int resp;

    if (sheet_modified) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 0);
	if (resp == GRETL_YES) {
	    get_data_from_sheet();
	} 
    }

    if (sheetcode == SHEET_NEW_DATASET) {
	empty_dataset_guard();
    }

    return FALSE;
}

void show_spreadsheet (SheetCode code) 
{
    GtkWidget *tmp, *button_box;
    GtkWidget *scroller, *main_vbox;
    GtkWidget *status_box, *mbar;
    GtkItemFactory *ifac;
    GtkAccelGroup *accel;

    if (gretlsheet != NULL) {
	gdk_window_raise(sheetwin->window);
	return;
    }

    if (datainfo->v == 1) {
	errbox(_("Please add a variable to the dataset first"));
	return;
    }

    sheetcode = code;
    sheetwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(sheetwin), _("gretl: edit data"));
    gtk_widget_set_usize(GTK_WIDGET(sheetwin), 600, 400);

    gtk_signal_connect(GTK_OBJECT(sheetwin), "delete_event",
		       (GtkSignalFunc) sheet_delete_event, 
		       NULL);

    main_vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 2); 
    gtk_container_add(GTK_CONTAINER(sheetwin), main_vbox);
    gtk_widget_show(main_vbox);

    /* add menu bar */
    accel = gtk_accel_group_new();
    ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", 
				accel);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(ifac, menu_translate, NULL, NULL);
#endif
    gtk_item_factory_create_items(ifac, 
				  sizeof sheet_items / sizeof sheet_items[0],
				  sheet_items, sheetwin);

    if (complex_subsampled()) {
	disable_obs_menu(ifac);
    }

    mbar = gtk_item_factory_get_widget(ifac, "<main>");
    gtk_accel_group_attach(accel, GTK_OBJECT (sheetwin));
    gtk_box_pack_start(GTK_BOX(main_vbox), mbar, FALSE, TRUE, 0);
    gtk_widget_show(mbar);

    status_box = gtk_hbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(status_box), 0);
    gtk_box_pack_start(GTK_BOX(main_vbox), status_box, FALSE, TRUE, 0);
    gtk_widget_show(status_box);

    locator = gtk_label_new(""); 
    gtk_widget_set_usize(locator, 160, 20);
    gtk_box_pack_start(GTK_BOX(status_box), locator, FALSE, TRUE, 0);
    gtk_widget_show(locator);

    topentry = gtk_entry_new(); 
    gtk_box_pack_start(GTK_BOX(status_box), topentry, TRUE, TRUE, 0); 
    gtk_widget_show(topentry);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(main_vbox), scroller, TRUE, TRUE, 0);
    gtk_widget_show(scroller);

    gretlsheet = gtk_sheet_new(datainfo->n, datainfo->v-1, "datasheet");
    gtk_container_add(GTK_CONTAINER(scroller), gretlsheet);
    gtk_widget_show(gretlsheet);

    /* apply and close buttons */
    button_box = gtk_hbox_new (FALSE, 5);
    gtk_box_set_homogeneous (GTK_BOX (button_box), TRUE);
    gtk_box_pack_start (GTK_BOX (main_vbox), button_box, FALSE, FALSE, 0);
    gtk_widget_show(button_box);

    gtk_signal_connect(GTK_OBJECT(sheetwin), "destroy",
		       GTK_SIGNAL_FUNC(free_spreadsheet), NULL);

    tmp = gtk_button_new_with_label(_("Apply Changes"));
    gtk_box_pack_start (GTK_BOX (button_box), tmp, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(get_data_from_sheet), NULL);
    gtk_widget_show(tmp);

    tmp = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start (GTK_BOX (button_box), tmp, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(maybe_exit_sheet), sheetwin);
    gtk_widget_show(tmp);

    gtk_signal_connect(GTK_OBJECT(gtk_sheet_get_entry(GTK_SHEET(gretlsheet))),
		       "changed", (GtkSignalFunc) show_entry, 
		       NULL);
    gtk_signal_connect(GTK_OBJECT(gretlsheet),
		       "activate", (GtkSignalFunc) activate_sheet_cell,
		       NULL);
    gtk_signal_connect(GTK_OBJECT(topentry),
		       "changed", (GtkSignalFunc) show_sheet_entry, 
		       NULL);
    gtk_signal_connect(GTK_OBJECT(topentry),
		       "activate", (GtkSignalFunc) activate_sheet_entry,
		       NULL);
    gtk_signal_connect(GTK_OBJECT(gretlsheet),
		       "button_press_event", (GtkSignalFunc) do_sheet_popup, 
		       NULL);
    gtk_signal_connect(GTK_OBJECT(gretlsheet),
		       "set_cell", (GtkSignalFunc) parse_numbers,
		       NULL); 

    GTK_SHEET_SET_FLAGS(gretlsheet, GTK_SHEET_AUTORESIZE);

    if (code == SHEET_EDIT_DATASET) {
	add_data_to_sheet(gretlsheet);
    } else {
	add_skel_to_sheet(gretlsheet);
    }

    gtk_widget_show(sheetwin);
    set_dataset_locked(TRUE);
}

/* mechanism for locking dataset against changes while editing */

static int locked;

static void set_dataset_locked (gboolean s)
{
    locked = s;
    flip(mdata->ifac, "/Sample", !s);
}

int dataset_locked (void)
{
    if (locked) {
	errbox(_("You can't do this while editing the dataset"));
    }

    return locked;
}
