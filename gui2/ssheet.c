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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

#include "gretl.h"
#include "ssheet.h"
#include <errno.h>
#include <ctype.h>

#define CELL_WIDTH 12
#define CELL_MAXLEN 23

typedef struct _spreadsheet spreadsheet;
typedef struct _sheet_cell sheet_cell;

struct _spreadsheet {
    GtkWidget *table;
    GtkWidget *win;
    GtkWidget *locator;
    GtkWidget *popup;
    GtkWidget **col_labels;
    GtkWidget **row_labels;
    sheet_cell ***cells;
    sheet_cell *active_cell;
    gchar location[20];
    gint cols, rows;
    gint n_scalars;
    guint cid;
    guint point;
};

struct _sheet_cell {
    int row;
    int col;
    double value;
    char text[CELL_MAXLEN + 1];
    GtkWidget *entry;
    spreadsheet *sheet;
};

enum {
    SHEET_AT_END,
    SHEET_AT_POINT
};

static void sheet_add_var_callback (gpointer data, guint code, GtkWidget *w);
static void sheet_add_obs_callback (gpointer data, guint where, GtkWidget *w);
static void get_data_from_sheet (GtkWidget *w, spreadsheet *sheet);
static gint get_data_col_width (void);
static void set_locator_label (spreadsheet *sheet, int row, int col);
#if 0
static void add_data_column (spreadsheet *sheet);
#endif

static GtkItemFactoryEntry sheet_items[] = {
    { N_("/_Observation"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Observation/_Append obs"), NULL, sheet_add_obs_callback, SHEET_AT_END, 
      NULL, GNULL },
    { N_("/Observation/_Insert obs"), NULL, sheet_add_obs_callback, SHEET_AT_POINT, 
      NULL, GNULL },
    { N_("/Add _Variable"), NULL, sheet_add_var_callback, 0, NULL, GNULL }
};

static int sheet_modified;

const char *sheet_cell_format = "%12g";

static int get_char_from_key (int key)
{
    if (key >= GDK_0 && key <= GDK_9) {
	return '0' + key - GDK_0;
    }
    
    if (key >= GDK_KP_0 && key <= GDK_KP_9) {
	return '0' + key - GDK_KP_0;
    }

    if (key == GDK_minus || key == GDK_KP_Subtract) {
	return '-';
    }

    /* add decimal point here */

    return '\0';
}

static int start_editing_cell (sheet_cell *cell, int keyval)
{
    gtk_editable_set_editable(GTK_EDITABLE(cell->entry), TRUE);

    if (gtk_editable_get_position(GTK_EDITABLE(cell->entry)) > 0) {
	return 0;
    } else {
	int c = get_char_from_key(keyval);
	char start[2];

	sprintf(start, "%c", c);
	gtk_entry_set_text(GTK_ENTRY(cell->entry), start);
	gtk_editable_set_position(GTK_EDITABLE(cell->entry), -1);
	return 1;
    }
}

static int finish_editing_cell (sheet_cell *cell)
{
    GtkWidget *w = cell->entry;
    const gchar *txt = gtk_entry_get_text(GTK_ENTRY(w));
    double val;

    val = atof(txt);

    if (val != cell->value) {
	sheet_modified = 1;
    }
    cell->value = val;
    sprintf(cell->text, sheet_cell_format, cell->value);
    gtk_entry_set_text(GTK_ENTRY(w), cell->text);
    gtk_editable_set_editable(GTK_EDITABLE(cell->entry), FALSE);

    return 0;
}

static void abort_editing_cell (sheet_cell *cell)
{
    sprintf(cell->text, sheet_cell_format, cell->value);
    gtk_entry_set_text(GTK_ENTRY(cell->entry), cell->text);
}

static void sheet_move_one_cell (sheet_cell *cell, int key)
{
    int i = cell->row;
    int j = cell->col;

    if (key == GDK_Down || key == GDK_Return) {
	if (++i < cell->sheet->rows) {
	    gtk_widget_grab_focus((cell->sheet->cells[j][i])->entry);
	}
    }

    if (key == GDK_Right) {
	if (++j < cell->sheet->cols) {
	    gtk_widget_grab_focus((cell->sheet->cells[j][i])->entry);
	}
    }

    if (key == GDK_Left) {
	if (--j >= 0) {
	    gtk_widget_grab_focus((cell->sheet->cells[j][i])->entry);
	}
    }

    if (key == GDK_Up) {
	if (--i >= 0) {
	    gtk_widget_grab_focus((cell->sheet->cells[j][i])->entry);
	}
    }	
}

static int deactivate_cell (sheet_cell *cell)
{
    if (cell != NULL) {
	if (gtk_editable_get_editable(GTK_EDITABLE(cell->entry))) {
	    if (finish_editing_cell(cell)) return 1;
	}
	gtk_entry_set_has_frame(GTK_ENTRY(cell->entry), FALSE);
	gtk_editable_set_editable(GTK_EDITABLE(cell->entry), FALSE);
    }
    return 0;
}

static void set_active_cell (sheet_cell *cell)
{
    if (cell != cell->sheet->active_cell) {
	if (deactivate_cell(cell->sheet->active_cell)) return; 
	cell->sheet->active_cell = cell;
	gtk_entry_set_has_frame(GTK_ENTRY(cell->entry), TRUE);
	gtk_editable_set_position(GTK_EDITABLE(cell->entry), 0);
	set_locator_label(cell->sheet, cell->row, cell->col);
    }
}

static gint cell_key_press (GtkWidget *widget,
			    GdkEventKey *key,
			    sheet_cell *cell)
{
    if (key->keyval == GDK_Return || key->keyval == GDK_Up ||
	key->keyval == GDK_Down) {
	int err = 0;

	if (gtk_editable_get_editable(GTK_EDITABLE(cell->entry))) {
	    err = finish_editing_cell(cell);
	}

	if (!err) {
	    sheet_move_one_cell(cell, key->keyval);
	}

	return TRUE;
    }

    if (gtk_editable_get_editable(GTK_EDITABLE(cell->entry))) {
	if (key->keyval == GDK_Escape) {
	    abort_editing_cell(cell);
	} 
	/* let the editing proceed */
	return FALSE;
    }

    if (key->keyval == GDK_Right || key->keyval == GDK_Left) {
	sheet_move_one_cell(cell, key->keyval);
	return TRUE;
    }

    return start_editing_cell(cell, key->keyval);
}

static void sheet_cell_set_value (spreadsheet *sheet, 
				  int row, int col, double val)
{
    sheet_cell *cell = sheet->cells[col][row];

    cell->value = val;
    if (cell->value == NADBL) {
	strcpy(cell->text, "NA");
	gtk_entry_set_text(GTK_ENTRY(cell->entry), "");
    } else {
	sprintf(cell->text, sheet_cell_format, cell->value);
	gtk_entry_set_text(GTK_ENTRY(cell->entry), cell->text);
    }
}

static sheet_cell *sheet_cell_new (spreadsheet *sheet, 
				   int row, int col)
{
    sheet_cell *cell;

    cell = malloc(sizeof *cell);
    if (cell == NULL) return NULL;

    cell->row = row;
    cell->col = col;

    cell->entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(cell->entry), CELL_MAXLEN);
    gtk_entry_set_has_frame(GTK_ENTRY(cell->entry), FALSE);
    gtk_editable_set_editable(GTK_EDITABLE(cell->entry), FALSE);
    gtk_entry_set_width_chars(GTK_ENTRY(cell->entry), CELL_WIDTH);

    cell->sheet = sheet;

    g_signal_connect(G_OBJECT(cell->entry), "key-press-event",
		     G_CALLBACK(cell_key_press), cell);

    return cell;
}

static void destroy_sheet (GtkWidget *widget, spreadsheet **psheet)
{
    spreadsheet *sheet = *psheet;
    int i, j;

    if (sheet->popup != NULL) {
	gtk_widget_destroy(sheet->popup);
    }

    if (sheet->cells != NULL) {
	for (j=0; j<sheet->cols; j++) {
	    if (sheet->cells[j] != NULL) {
		for (i=0; i<sheet->rows; i++) {
		    free(sheet->cells[j][i]);
		}
		free(sheet->cells[j]);
	    }
	}
	free(sheet->cells);
	sheet->cells = NULL;
    }

    free(sheet);
    *psheet = NULL;
}

static gint click_cell (GtkWidget *entry,
			GdkEventButton *event,
			sheet_cell *cell)
{
    if (cell != cell->sheet->active_cell) {
	set_active_cell(cell);
	gtk_widget_grab_focus(cell->entry);
	return TRUE;
    }
    return FALSE;
}

static gint focus_in_cell (GtkWidget *widget,
			   GdkEventFocus *event,
			   sheet_cell *cell)
{
    set_active_cell(cell);
    
    return 0;
}

static spreadsheet *sheet_new (int rows, int cols, int nscalars)
{
    spreadsheet *sheet;
    int i, j;

    sheet = malloc(sizeof *sheet);
    if (sheet == NULL) return NULL;

    sheet->rows = rows;
    sheet->cols = cols;
    sheet->n_scalars = nscalars;

    sheet->col_labels = NULL;
    sheet->row_labels = NULL;

    sheet->cells = malloc(cols * sizeof *sheet->cells);
    if (sheet->cells == NULL) {
	free(sheet);
	return NULL;
    }

    for (j=0; j<cols; j++) {
	sheet->cells[j] = malloc(rows * sizeof **sheet->cells);
	if (sheet->cells[j] == NULL) {
	    for (i=0; i<j; i++) {
		free(sheet->cells[i]);
	    }
	    free(sheet->cells);
	    free(sheet);
	    return NULL;
	} 
	for (i=0; i<rows; i++) {
	    sheet->cells[j][i] = sheet_cell_new(sheet, i, j); 
	    if (sheet->cells[j][i] == NULL) {
		;
	    }
	}
    }

    sheet->col_labels = malloc(cols * sizeof *sheet->col_labels);
    if (sheet->col_labels == NULL) {
	destroy_sheet(NULL, &sheet);
	return NULL;
    }
    for (j=0; j<cols; j++) {
	sheet->col_labels[j] = gtk_label_new("");
    }

    sheet->row_labels = malloc(rows * sizeof *sheet->row_labels);
    if (sheet->row_labels == NULL) {
	destroy_sheet(NULL, &sheet);
	return NULL;
    }
    for (i=0; i<rows; i++) {
	sheet->row_labels[i] = gtk_label_new("");
    }    

    sheet->table = gtk_table_new(rows + 1, cols + 1, TRUE);
    if (sheet->table == NULL) {
	destroy_sheet(NULL, &sheet);
	return NULL;
    }  

    /* attach all cells and labels to table */
    for (i=0; i<rows; i++) {
	gtk_table_attach(GTK_TABLE(sheet->table), 
			 sheet->row_labels[i], 
			 0, 1, i+1, i+2,
			 0, 0, 0, 0);
	for (j=0; j<cols; j++) {
	    gtk_table_attach(GTK_TABLE(sheet->table), 
			     (sheet->cells[j][i])->entry, 
			     j+1, j+2, i+1, i+2,
			     0, 0, 0, 0);
	    g_signal_connect(G_OBJECT((sheet->cells[j][i])->entry),
			     "button-press-event", 
			     G_CALLBACK(click_cell), sheet->cells[j][i]);
	    g_signal_connect(G_OBJECT((sheet->cells[j][i])->entry),
			     "focus-in-event", 
			     G_CALLBACK(focus_in_cell), sheet->cells[j][i]);

	}
    } 

    for (j=0; j<cols; j++) {
	gtk_table_attach(GTK_TABLE(sheet->table), 
			 sheet->col_labels[j], 
			 j+1, j+2, 0, 1,
			 0, 0, 0, 0);
    }

    sheet->win = NULL;
    sheet->locator = NULL;
    sheet->popup = NULL;
    sheet->active_cell = NULL;
    sheet->cid = 0;

    return sheet;
}

static void label_sheet_column (spreadsheet *sheet, int col,
				const char *str)
{
    gtk_label_set_text(GTK_LABEL(sheet->col_labels[col]), str);
}

static void label_sheet_row (spreadsheet *sheet, int row,
			     const char *str)
{
    gtk_label_set_text(GTK_LABEL(sheet->row_labels[row]), str);
}

const gchar *sheet_get_column_label (spreadsheet *sheet, int col)
{
    return gtk_label_get_text(GTK_LABEL(sheet->col_labels[col]));
}

const gchar *sheet_get_row_label (spreadsheet *sheet, int row)
{
    return gtk_label_get_text(GTK_LABEL(sheet->row_labels[row]));
}

/* .................................................................. */

static void set_locator_label (spreadsheet *sheet, int row, int col)
{
    const char *col_label;
    const char *row_label;

    col_label = gtk_label_get_text(GTK_LABEL(sheet->col_labels[col]));
    row_label = gtk_label_get_text(GTK_LABEL(sheet->row_labels[row]));
    
    sprintf(sheet->location, "%s, %s", col_label, row_label);
    gtk_statusbar_pop(GTK_STATUSBAR(sheet->locator), sheet->cid);
    gtk_statusbar_push(GTK_STATUSBAR(sheet->locator), 
		       sheet->cid, sheet->location);
}

static void real_add_new_var (spreadsheet *sheet, const char *varname)
{
    return;
}

static void real_add_new_obs (spreadsheet *sheet, const char *marker)
{
    return;
}

/* ........................................................... */

static void name_new_var (GtkWidget *widget, dialog_t *ddata) 
{
    spreadsheet *sheet = (spreadsheet *) ddata->data;
    const gchar *buf;
    char varname[9];

    buf = gtk_entry_get_text (GTK_ENTRY (ddata->edit));

    if (blank_entry(buf, ddata)) return;
    if (validate_varname(buf)) return;

    *varname = 0;
    strncat(varname, buf, 8);

    close_dialog(ddata);
    real_add_new_var(sheet, varname);
}

/* ........................................................... */

static void name_new_obs (GtkWidget *widget, dialog_t *ddata) 
{
    spreadsheet *sheet = (spreadsheet *) ddata->data;
    const gchar *buf;
    char obsmarker[9];

    buf = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (blank_entry(buf, ddata)) return;

    *obsmarker = 0;
    strncat(obsmarker, buf, 8);

    close_dialog(ddata);
    real_add_new_obs(sheet, obsmarker);
}

/* ........................................................... */

static void name_var_dialog (spreadsheet *sheet) 
{
    edit_dialog (_("gretl: name variable"), 
		 _("Enter name for new variable\n"
		   "(max. 8 characters)"),
		 NULL, name_new_var, sheet, 
		 0, 0);
}

/* ........................................................... */

static void new_case_dialog (spreadsheet *sheet) 
{
    edit_dialog (_("gretl: case marker"), 
		 _("Enter case marker for new obs\n"
		   "(max. 8 characters)"),
		 NULL, name_new_obs, sheet, 
		 0, 0);
}

/* ........................................................... */

static void sheet_add_var_callback (gpointer data, guint u, GtkWidget *w)
{
    spreadsheet *sheet = (spreadsheet *) data;

    name_var_dialog(sheet);
}

/* ........................................................... */

static void sheet_add_obs_callback (gpointer data, guint where, GtkWidget *w)
{
    spreadsheet *sheet = (spreadsheet *) data;

    sheet->point = where;

    if (datainfo->markers) {
	new_case_dialog(sheet);
    } else {
	real_add_new_obs(sheet, "");
    }
}

/* ........................................................... */

static void popup_sheet_add_obs (GtkWidget *w, spreadsheet *sheet)
{
    sheet_add_obs_callback(sheet, SHEET_AT_END, NULL);
}

static void popup_sheet_insert_obs (GtkWidget *w, spreadsheet *sheet)
{
    sheet_add_obs_callback(sheet, SHEET_AT_POINT, NULL);
}

static void popup_sheet_add_var (GtkWidget *w, spreadsheet *sheet)
{
    sheet_add_var_callback(sheet, 0, NULL);
}

/* ........................................................... */

static void build_sheet_popup (GtkWidget **popup, spreadsheet *sheet)
{
    if (*popup != NULL) return;

    *popup = gtk_menu_new();

    add_popup_item(_("Add Variable"), *popup, 
		   G_CALLBACK(popup_sheet_add_var),
		   sheet);
    add_popup_item(_("Add Observation"), *popup,
		   G_CALLBACK(popup_sheet_add_obs),
		   sheet);
    add_popup_item(_("Insert Observation"), *popup,
		   G_CALLBACK(popup_sheet_insert_obs),
		   sheet);
}

/* ........................................................... */

static void get_data_from_sheet (GtkWidget *w, spreadsheet *sheet)
{
    gint i, t, n = datainfo->n, oldv = datainfo->v; 
    gint orig_cols, newvars, newobs, missobs = 0;
    gint colnum;

    newobs = sheet->rows - n;
    orig_cols = datainfo->v - 1 - sheet->n_scalars;
    newvars = sheet->cols - orig_cols;

    if (newobs > 0) {
	if (grow_nobs(newobs, &Z, datainfo)) {
	    errbox(_("Failed to allocate memory for new data"));
	    return;
	}
	n = datainfo->n;
    }

    if (newvars > 0) {
	if (dataset_add_vars(newvars, &Z, datainfo)) {
	    errbox(_("Failed to allocate memory for new data"));
	    return;
	}
	for (i=0; i<newvars; i++) { 
	    const gchar *newname;

	    newname = sheet_get_column_label(sheet, orig_cols + 1 + i); /* FIXME? */
	    strcpy(datainfo->varname[i + oldv], newname);
	    strcpy(VARLABEL(datainfo, i + oldv), "");
	}
    }

    colnum = 0;
    for (i=1; i<datainfo->v; i++) {
	if (datainfo->vector[i] == 0) continue;
	colnum++;
	for (t=0; t<n; t++) {
	    Z[i][t] = (sheet->cells[colnum][t])->value;
	}
    }

    if (datainfo->markers && datainfo->S != NULL) {
	for (t=0; t<n; t++) {
	    strcpy(datainfo->S[t], sheet_get_row_label(sheet, t));
	}
    }

    data_status |= (GUI_DATA | MODIFIED_DATA);
    register_data(NULL, NULL, 0);

    if (missobs) {
	infobox(_("Warning: there were missing observations"));
    } else {
	infobox(_("Data updated OK"));
    }

    sheet_modified = 0;
}

/* ........................................................... */

static void add_data_to_sheet (spreadsheet *sheet)
{
    char rowlabel[9];
    int i, t, colnum;

    /* insert the variable names */
    colnum = 0;
    for (i=1; i<datainfo->v; i++) {
	if (datainfo->vector[i] == 0) continue;
	label_sheet_column(sheet, colnum, datainfo->varname[i]);
	colnum++;
    }

    /* insert observation marker column */
    for (t=0; t<sheet->rows; t++) {
	if (datainfo->markers) {
	    strcpy(rowlabel, datainfo->S[t]);
	} else {
	    ntodate(rowlabel, t, datainfo);
	}
	label_sheet_row(sheet, t, rowlabel);
    }

    /* insert the data values */
    for (t=0; t<sheet->rows; t++) {
	colnum = 0;
	for (i=1; i<datainfo->v; i++) {
	    /* don't put scalars into the spreadsheet */
	    if (datainfo->vector[i] == 0) continue;
	    sheet_cell_set_value(sheet, t, colnum, Z[i][t]);
	    colnum++;
	}
    }
}

/* ........................................................... */

static void add_skel_to_sheet (spreadsheet *sheet)
{
    char rowlabel[9];
    int i, t;

    /* insert the variable names */
    for (i=1; i<datainfo->v; i++) {
	label_sheet_column(sheet, i - 1, datainfo->varname[i]);
    }

    /* insert observation markers */
    for (t=0; t<sheet->rows; t++) {
	ntodate(rowlabel, t, datainfo);
	label_sheet_row(sheet, t, rowlabel);
    }

    /* insert "missing" data values */
    for (t=0; t<sheet->rows; t++) {
	for (i=0; i<sheet->cols; i++) {
	    sheet_cell_set_value(sheet, t, i, NADBL);
	}
    }
}

/* ........................................................... */

static gint get_string_width (const gchar *str)
{
    gint width;

    GtkWidget *w;
    PangoLayout *pl;
    PangoContext *pc;

    w = gtk_label_new(NULL);
    pc = gtk_widget_get_pango_context(w);

    pl = pango_layout_new(pc);
    pango_layout_set_text(pl, str, -1);
    pango_layout_get_pixel_size(pl, &width, NULL);

    gtk_widget_destroy(w);
    g_object_unref(G_OBJECT(pl));

    return width;
}

static gint get_obs_col_width (void)
{
    static gint width;

    if (width == 0) width = get_string_width("XXXXXXXXX");
    return width;
}

static gint get_data_col_width (void)
{
    static gint width;

    if (width == 0) width = get_string_width("-000000.000000");
    return width;
}

/* ........................................................... */

static gint maybe_exit_sheet (GtkWidget *w, spreadsheet *sheet)
{
    int resp;

    if (sheet_modified) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    get_data_from_sheet(NULL, sheet);
	}
	else if (resp == GRETL_CANCEL || resp == -1) return FALSE;
    }
  
    gtk_widget_destroy(sheet->win);

    return FALSE;
}

/* ........................................................... */

static int get_n_scalars (const DATAINFO *pdinfo)
{
    int i, ns = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->vector[i] == 0) ns++;
    }

    return ns;
}

/* ........................................................... */

void show_spreadsheet (DATAINFO *pdinfo) 
{
    static spreadsheet *sheet;    
    GtkWidget *tmp, *button_box;
    GtkWidget *scroller, *main_vbox;
    GtkWidget *status_box, *mbar;
    GtkItemFactory *ifac;
    int sheetwidth;
    int nscalars;

    if (sheet != NULL) {
	gdk_window_raise(sheet->win->window);
	return;
    }

    if (pdinfo == NULL && datainfo->v == 1) {
	errbox(_("Please add a variable to the dataset first"));
	return;
    }

    if (pdinfo == NULL) {
	nscalars = get_n_scalars(datainfo);
	sheet = sheet_new(datainfo->n, datainfo->v - 1 - nscalars, nscalars);
    } else {
	sheet = sheet_new(pdinfo->n, 1, 0);
    }

    if (sheet == NULL) return;   

    sheetwidth = get_obs_col_width() + 6 * get_data_col_width() + 14;

    sheet->win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(sheet->win), _("gretl: edit data"));
    gtk_window_set_default_size(GTK_WINDOW(sheet->win), sheetwidth, 400);

    main_vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 5); 
    gtk_container_add(GTK_CONTAINER(sheet->win), main_vbox);

    /* add menu bar */
    ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", 
				NULL);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(ifac, menu_translate, NULL, NULL);
#endif
    gtk_item_factory_create_items(ifac, 
				  sizeof sheet_items / sizeof sheet_items[0],
				  sheet_items, sheet);
    mbar = gtk_item_factory_get_widget(ifac, "<main>");
    gtk_box_pack_start(GTK_BOX(main_vbox), mbar, FALSE, FALSE, 0);

    build_sheet_popup(&sheet->popup, sheet);

    status_box = gtk_hbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(status_box), 0);
    gtk_box_pack_start(GTK_BOX(main_vbox), status_box, FALSE, FALSE, 0);

    sheet->locator = gtk_statusbar_new(); 
    gtk_widget_set_size_request(sheet->locator, 2 * get_obs_col_width(), 20);
    gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(sheet->locator), FALSE);
    sheet->cid = gtk_statusbar_get_context_id (GTK_STATUSBAR(sheet->locator), 
					       "current row and column");
    gtk_box_pack_start(GTK_BOX(status_box), sheet->locator, FALSE, FALSE, 0);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
                                    GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW(scroller), 
					 GTK_SHADOW_NONE);

    gtk_box_pack_start(GTK_BOX(main_vbox), scroller, TRUE, TRUE, 0);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller),
					  sheet->table);

    /* apply and close buttons */
    button_box = gtk_hbox_new (FALSE, 5);
    gtk_box_set_homogeneous (GTK_BOX (button_box), TRUE);
    gtk_box_pack_start (GTK_BOX(main_vbox), button_box, FALSE, FALSE, 0);

    g_signal_connect(G_OBJECT(sheet->win), "destroy",
		     G_CALLBACK(destroy_sheet), &sheet);

    tmp = gtk_button_new_from_stock(GTK_STOCK_APPLY);
    gtk_box_pack_start (GTK_BOX (button_box), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(get_data_from_sheet), sheet);

    tmp = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_box_pack_start (GTK_BOX (button_box), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(maybe_exit_sheet), sheet);

    if (pdinfo != NULL) {
	add_skel_to_sheet(sheet);
    } else {
	add_data_to_sheet(sheet);
    }

    sheet_modified = 0;

    gtk_widget_show_all(sheet->win);
    set_active_cell(sheet->cells[0][0]);
}

