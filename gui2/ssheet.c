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
#include "treeutils.h"
#include "ssheet.h"
#include <errno.h>
#include <ctype.h>

#define SHEET_PRECISION 10

typedef struct {
    GtkWidget *view;
    GtkWidget *win;
    GtkWidget *locator;
    GtkWidget *popup;
    gint datacols, datarows;
    gint padcols, totcols;
    gint n_scalars;
    guint cid;
    guint point;
} spreadsheet;

enum {
    SHEET_AT_END,
    SHEET_AT_POINT
};

static void sheet_add_var_callback (gpointer data, guint code, GtkWidget *w);
static void sheet_add_obs_callback (gpointer data, guint where, GtkWidget *w);
static void get_data_from_sheet (GtkWidget *w, spreadsheet *sheet);
static GtkCellRenderer *datacell_new (GtkListStore *store, gint colnum);
static void set_up_sheet_column (GtkTreeViewColumn *column, gint width);
static gint get_data_col_width (void);
static void add_data_column (spreadsheet *sheet);

static GtkItemFactoryEntry sheet_items[] = {
    { N_("/_Observation"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Observation/_Append obs"), NULL, sheet_add_obs_callback, SHEET_AT_END, 
      NULL, GNULL },
    { N_("/Observation/_Insert obs"), NULL, sheet_add_obs_callback, SHEET_AT_POINT, 
      NULL, GNULL },
    { N_("/Add _Variable"), NULL, sheet_add_var_callback, 0, NULL, GNULL }
};

/* .................................................................. */

static gint sheet_cell_edited (GtkCellRendererText *cell,
			       const gchar *path_string,
			       const gchar *new_text,
			       gpointer data)
{
    GtkTreeModel *model = (GtkTreeModel *) data;
    gint err = 0;

    err = check_atof(new_text);

    if (err) {
	errbox(get_gretl_errmsg());
    } else {
	gint *column;
	GtkTreePath *path;
	GtkTreeIter iter;

	column = g_object_get_data(G_OBJECT(cell), "column");
	path = gtk_tree_path_new_from_string (path_string);
	gtk_tree_model_get_iter(model, &iter, path);
	gtk_list_store_set(GTK_LIST_STORE(model), &iter, 
			   GPOINTER_TO_INT(column), new_text, -1);
	gtk_tree_path_free(path);
    }

    return FALSE;
}

/* .................................................................. */

static void real_add_new_var (spreadsheet *sheet, const char *varname)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view); 
    GtkTreeModel *model;
    GtkTreeViewColumn *column;
    GtkCellRenderer *renderer;
    GList *collist = NULL;
    gint i, oldcols, cols;

    oldcols = sheet->totcols;

    add_data_column(sheet);

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 1, NULL);
    g_object_set(renderer, "xalign", 1.0, NULL);
	
    if (sheet->totcols == oldcols) {
	column = gtk_tree_view_get_column(view, sheet->datacols);
	gtk_tree_view_remove_column(view, column);
    }

    column = gtk_tree_view_column_new();
    gtk_tree_view_column_set_title(column, varname);
    cols = gtk_tree_view_insert_column(view, column, sheet->datacols);
	
    set_up_sheet_column(column, get_data_col_width()); 

    model = gtk_tree_view_get_model(view);

    collist = gtk_tree_view_get_columns(view);
    for (i=0; i<sheet->totcols-2; i++) {
	column = GTK_TREE_VIEW_COLUMN(collist->data);
	gtk_tree_view_column_clear(column);

	if (i > 0 && i <= sheet->datacols) {
	    GtkCellRenderer *datacell;

	    datacell = datacell_new(GTK_LIST_STORE(model), i);
	    gtk_tree_view_column_pack_start(column, datacell, TRUE);
	    gtk_tree_view_column_set_attributes (column,
						 datacell,
						 "text", i, 
						 "editable", sheet->totcols - 1, 
						 NULL);
	} else { /* non-data cells */
	    gtk_tree_view_column_pack_start(column, renderer, TRUE);
	    gtk_tree_view_column_set_attributes (column,
						 renderer,
						 "text", i, 
						 "editable", sheet->totcols - 2,
						 NULL);
	}
	collist = collist->next;
    }
    if (collist) g_list_free(collist);
}

/* .................................................................. */

static void real_add_new_obs (spreadsheet *sheet, const char *obsname)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    gint pointpath = 0;
    GtkListStore *store;
    GtkTreeIter iter;
    gchar rowlabel[10];
    gint i;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));

    if (sheet->point == SHEET_AT_END) {
	gtk_list_store_append(store, &iter);
    } 
    else if (sheet->point == SHEET_AT_POINT) {
	GtkTreePath *path;
	GtkTreeViewColumn *column;

	gtk_tree_view_get_cursor(view, &path, &column);
	gtk_tree_model_get_iter(GTK_TREE_MODEL(store), &iter, path);
	pointpath = gtk_tree_path_get_indices(path)[0];
	gtk_list_store_insert(store, &iter, pointpath);
	gtk_tree_path_free(path);
    } 
    else return;

    sheet->datarows += 1;

    if (datainfo->markers) {
	gtk_list_store_set(store, &iter, 0, obsname, -1);
    } else if (sheet->point == SHEET_AT_END) {
	ntodate(rowlabel, sheet->datarows - 1, datainfo);
	gtk_list_store_set(store, &iter, 0, rowlabel, -1);
    }

    for (i=1; i<=sheet->datacols; i++) {
	gtk_list_store_set(store, &iter, i, "", -1);
    }

    for (i=1; i<=sheet->padcols; i++) {
	gtk_list_store_set(store, &iter, sheet->datacols + i, "", -1);
    }

    gtk_list_store_set(store, &iter, 
		       sheet->totcols - 2, FALSE, 
		       sheet->totcols - 1, TRUE,
		       -1);

    if (sheet->point == SHEET_AT_POINT && !datainfo->markers) {
	for (i=pointpath; i<sheet->datarows; i++) {
	    ntodate(rowlabel, i, datainfo);
	    gtk_list_store_set(store, &iter, 0, rowlabel, -1);
	    gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	}
    } 

    if (sheet->point == SHEET_AT_END) {
	GtkTreePath *path;
	GtkTreeViewColumn *column;
	GtkAdjustment *adj;
	gchar *pathstr;
	gdouble adjval;

	pathstr = g_strdup_printf("%d", sheet->datarows - 1);
	path = gtk_tree_path_new_from_string(pathstr);
	column = gtk_tree_view_get_column(view, 1);
	gtk_tree_view_set_cursor(view, path, column, TRUE);
	adj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(gtk_widget_get_ancestor
								      (GTK_WIDGET(view),
								       GTK_TYPE_BIN)));
	adjval = gtk_adjustment_get_value(adj);
	gtk_adjustment_set_value (adj, adjval + 30); /* why is this hack needed? */
	gtk_tree_path_free(path);
	g_free(pathstr);
    }
}

/* ........................................................... */

static void name_new_var (GtkWidget *widget, dialog_t *ddata) 
{
    spreadsheet *sheet = (spreadsheet *) ddata->data;
    const gchar *buf;
    char newvarname[9];

    buf = gtk_entry_get_text (GTK_ENTRY (ddata->edit));

    if (blank_entry(buf, ddata)) return;
    if (validate_varname(buf)) return;

    *newvarname = 0;
    strncat(newvarname, buf, 8);

    close_dialog(ddata);
    real_add_new_var(sheet, newvarname);
}

/* ........................................................... */

static void name_new_obs (GtkWidget *widget, dialog_t *ddata) 
{
    spreadsheet *sheet = (spreadsheet *) ddata->data;
    const gchar *buf;
    char newobsmarker[9];

    buf = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (blank_entry(buf, ddata)) return;

    *newobsmarker = 0;
    strncat(newobsmarker, buf, 8);

    close_dialog(ddata);
    real_add_new_obs(sheet, newobsmarker);
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

static void add_data_column (spreadsheet *sheet)
{
    GType *types;
    GtkListStore *old_store, *new_store;
    GtkTreeIter old_iter, new_iter;
    gint i, row, newcolnum;

    sheet->datacols += 1;

    if (sheet->padcols > 0) sheet->padcols -= 1;
    else sheet->totcols += 1;

    /* make an expanded column types list */
    types = mymalloc(sheet->totcols * sizeof *types);

    /* configure the types */
    types[0] = G_TYPE_STRING;
    for (i=1; i<=sheet->datacols; i++) types[i] = G_TYPE_STRING;
    for (i=0; i<sheet->padcols; i++) types[sheet->datacols+1+i] = G_TYPE_STRING;
    for (i=sheet->totcols-2; i<=sheet->totcols-1; i++) types[i] = G_TYPE_BOOLEAN;

    newcolnum = sheet->datacols;

    old_store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sheet->view)));
    new_store = gtk_list_store_newv (sheet->totcols, types);
    free(types);

    /* go to start of old and new lists */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(old_store), &old_iter);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(new_store), &new_iter);

    /* copy across the old data */
    for (row=0; row<sheet->datarows; row++) {
	gchar *str;

	gtk_list_store_append(new_store, &new_iter);

	/* column 0: markers */
	gtk_tree_model_get(GTK_TREE_MODEL(old_store), &old_iter, 0, &str, -1);
	gtk_list_store_set(new_store, &new_iter, 0, str, -1);
	g_free(str);

	/* original data values */
	for (i=1; i<newcolnum; i++) {
	    gtk_tree_model_get(GTK_TREE_MODEL(old_store), &old_iter, i, &str, -1);
	    gtk_list_store_set(new_store, &new_iter, i, str, -1);
	    g_free(str);
	}

	/* new data values blank */
	gtk_list_store_set(new_store, &new_iter, newcolnum, "", -1);

	/* any padding cols */
	for (i=1; i<=sheet->padcols; i++) {
	    gtk_list_store_set(new_store, &new_iter, newcolnum + i, "", -1);
	}	

	/* editable flags */
	gtk_list_store_set(new_store, &new_iter, 
			   sheet->totcols - 2, FALSE, 
			   sheet->totcols - 1, TRUE, -1);

	gtk_tree_model_iter_next(GTK_TREE_MODEL(old_store), &old_iter);
    }    

    gtk_tree_view_set_model(GTK_TREE_VIEW(sheet->view), GTK_TREE_MODEL(new_store));
    g_object_unref(G_OBJECT(new_store));
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

/* ......................................................... */

static gboolean update_selected (GtkTreeSelection *selection, spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreePath *path;
    GtkTreeViewColumn *column;

    gtk_tree_view_get_cursor(view, &path, &column);

    if (path && column) {
	GtkTreeModel *model;
	GtkTreeIter iter;
	gchar *row_label;
	    
	model = gtk_tree_view_get_model(view);
	gtk_tree_model_get_iter(model, &iter, path);
	gtk_tree_model_get(model, &iter, 0, &row_label, -1);
	gtk_statusbar_pop(GTK_STATUSBAR(sheet->locator), sheet->cid);
	gtk_statusbar_push(GTK_STATUSBAR(sheet->locator), 
			   sheet->cid, row_label);
	g_free(row_label);
	gtk_tree_path_free(path);
    }

    return FALSE;
}

/* ........................................................... */

static void get_data_from_sheet (GtkWidget *w, spreadsheet *sheet)
{
    gint i, t, n = datainfo->n, oldv = datainfo->v; 
    gint orig_cols, newvars, newobs, missobs = 0;
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkTreeViewColumn *column;
    GtkTreeModel *model;
    gint colnum;

    newobs = sheet->datarows - n;
    orig_cols = datainfo->v - 1 - sheet->n_scalars;
    newvars = sheet->datacols - orig_cols;

    model = gtk_tree_view_get_model(view);

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

	    column = gtk_tree_view_get_column(view, orig_cols + 1 + i);
	    newname = gtk_tree_view_column_get_title(column);
	    strcpy(datainfo->varname[i + oldv], newname);
	    strcpy(VARLABEL(datainfo, i + oldv), "");
	}
    }

    colnum = 0;
    for (i=1; i<datainfo->v; i++) {
	if (datainfo->vector[i] == 0) continue;
	colnum++;
	gtk_tree_model_get_iter_first(model, &iter);	
	for (t=0; t<n; t++) {
	    gchar *numstr;

	    gtk_tree_model_get(model, &iter, colnum, &numstr, -1);
	    if (*numstr) {
		Z[i][t] = atof(numstr); 
	    } else {
		Z[i][t] = NADBL;
		missobs = 1;
	    }
	    g_free(numstr);
	    gtk_tree_model_iter_next(model, &iter);
	}
    }

    if (datainfo->markers && datainfo->S != NULL) {
	gtk_tree_model_get_iter_first(model, &iter);
	for (t=0; t<n; t++) {
	    gchar *marker;

	    gtk_tree_model_get(model, &iter, 0, &marker, -1);
	    strcpy(datainfo->S[t], marker);
	    g_free(marker);
	    gtk_tree_model_iter_next(model, &iter);
	}
    }

    data_status |= (GUI_DATA|MODIFIED_DATA);
    register_data(NULL, NULL, 0);

    if (missobs)
	infobox(_("Warning: there were missing observations"));
    else
	infobox(_("Data updated OK"));
}

/* ........................................................... */

static void add_data_to_sheet (spreadsheet *sheet)
{
    gchar rowlabel[10];
    gint i, t, colnum, n = datainfo->n;
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkListStore *store;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));

    /* insert observation markers */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	if (datainfo->markers) 
	    strcpy(rowlabel, datainfo->S[t]);
	else
	    ntodate(rowlabel, t, datainfo);
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, rowlabel, -1);
    }

    sheet->datarows = t;

    /* insert the data values */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	char numstr[32];

	colnum = 0;
	for (i=1; i<datainfo->v; i++) {
	    /* don't put scalars into the spreadsheet */
	    if (datainfo->vector[i] == 0) continue;
	    if (na(Z[i][t])) {
		*numstr = '\0';
	    } else {
		sprintf(numstr, "%.*g", SHEET_PRECISION, Z[i][t]);
	    }
	    colnum++;
	    gtk_list_store_set(store, &iter, colnum, numstr, -1);
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }

    /* insert padding cols if needed */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	for (i=0; i<sheet->padcols; i++) {
	    gtk_list_store_set(store, &iter, i + datainfo->v, "", -1);
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }

    /* set the editable flags */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	gtk_list_store_set(store, &iter, 
			   sheet->totcols - 2, FALSE,
			   sheet->totcols - 1, TRUE, -1);
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }
}

/* ........................................................... */

static void add_skel_to_sheet (spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkListStore *store;
    gchar rowlabel[10];
    gint i, t, n = datainfo->n;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));

    /* insert observation markers */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	ntodate(rowlabel, t, datainfo);
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, rowlabel, -1);
    }

    sheet->datarows = t;

    /* insert "missing" data values */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	for (i=1; i<datainfo->v; i++) {
	    gtk_list_store_set(store, &iter, i, "", -1);
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }

    /* insert padding cols if needed */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	for (i=0; i<sheet->padcols; i++) {
	    gtk_list_store_set(store, &iter, i + datainfo->v, "", -1);
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }

    /* set the editable flags */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=0; t<n; t++) {
	gtk_list_store_set(store, &iter, 
			   sheet->totcols - 2, FALSE,
			   sheet->totcols - 1, TRUE, -1);
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
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

static void set_up_sheet_column (GtkTreeViewColumn *column, gint width)
{
    gtk_tree_view_column_set_alignment(column, 0.5); /* header centered */
    gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(column, width);
}

/* ........................................................... */

static GtkCellRenderer *datacell_new (GtkListStore *store, gint colnum)
{
    GtkCellRenderer *datacell;

    datacell = gtk_cell_renderer_text_new();
    g_object_set(datacell, "ypad", 1, NULL);
    g_object_set(datacell, "xalign", 1.0, NULL);
    g_signal_connect (G_OBJECT (datacell), "edited",
		      G_CALLBACK (sheet_cell_edited), GTK_TREE_MODEL(store));
    g_object_set_data(G_OBJECT(datacell), "column", (gint *) colnum);

    return datacell;
}

/* ........................................................... */

static GtkWidget *data_sheet_new (spreadsheet *sheet, gint nobs, gint nvars)
{
    GtkListStore *store; 
    GtkWidget *view;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    GType *types;
    gint i, width, colnum;

    sheet->datacols = nvars - 1; /* don't show the constant */

    /* we'll drop any scalar variables from the spreadsheet */
    for (i=1; i<datainfo->v; i++) {
	if (datainfo->vector[i] == 0) {
	    sheet->datacols -= 1;
	    sheet->n_scalars += 1;
	}
    }

    if (sheet->datacols < 6) sheet->padcols = 6 - sheet->datacols;
    else sheet->padcols = 0;

    /* obs, data, padding, boolean cols */
    sheet->totcols = 1 + sheet->datacols + sheet->padcols + 2;

    types = mymalloc(sheet->totcols * sizeof *types);
    if (types == NULL) return NULL;

    types[0] = G_TYPE_STRING;                             /* observation marker */
    for (i=1; i<=sheet->datacols; i++) 
	types[i] = G_TYPE_STRING;                         /* string rep. of data values */
    for (i=0; i<sheet->padcols; i++) 
	types[i + sheet->datacols + 1] = G_TYPE_STRING;   /* padding columns */
    types[sheet->totcols - 2] = G_TYPE_BOOLEAN;           /* FALSE editable flag */
    types[sheet->totcols - 1] = G_TYPE_BOOLEAN;           /* TRUE editable flag */

    store = gtk_list_store_newv (sheet->totcols, types);
    free(types);

    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    gtk_tree_view_set_rules_hint (GTK_TREE_VIEW(view), TRUE);

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 1, NULL);
    g_object_set(renderer, "xalign", 1.0, NULL);

    /* construct the observation marker column */
    width = get_obs_col_width();
    column = gtk_tree_view_column_new_with_attributes (NULL,
						       renderer,
						       "text", 0, 
						       "editable", 
						       sheet->totcols - 2,
						       NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
    set_up_sheet_column(column, width);

    /* construct the data columns */
    width = get_data_col_width();
    colnum = 0;
    for (i=1; i<nvars; i++) {
	GtkCellRenderer *datacell;

	if (datainfo->vector[i] == 0) continue;

	colnum++;
	datacell = datacell_new(store, colnum);
	column = gtk_tree_view_column_new_with_attributes (datainfo->varname[i],
							   datacell,
							   "text", 
							   colnum, 
							   "editable", 
							   sheet->totcols - 1,
							   NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	set_up_sheet_column(column, width);
    }

    /* add some padding columns if needed */
    for (i=0; i<sheet->padcols; i++) {
	column = gtk_tree_view_column_new_with_attributes (NULL,
							   renderer,
							   "text", 
							   i + sheet->datacols + 1,
							   "editable", 
							   sheet->totcols - 2,
							   NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	set_up_sheet_column(column, width);
    }

    /* set the selection properties on the tree view */
    select = gtk_tree_view_get_selection (GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);

    g_signal_connect (G_OBJECT(select), "changed",
		      G_CALLBACK(update_selected), sheet);

    return view;
}

/* ........................................................... */

static void free_spreadsheet (GtkWidget *widget, spreadsheet **psheet) 
{
    spreadsheet *sheet = *psheet;

    gtk_widget_destroy(sheet->popup);
    free(sheet);
    *psheet = NULL;
}

/* ........................................................... */

static spreadsheet *sheet_new (void)
{
    spreadsheet *sheet;

    sheet = malloc(sizeof *sheet);
    if (sheet == NULL) return NULL;

    sheet->view = NULL;
    sheet->win = NULL;
    sheet->locator = NULL;
    sheet->popup = NULL;
    sheet->datacols = sheet->datarows = 0;
    sheet->padcols = sheet->totcols = 0;
    sheet->n_scalars = 0;
    sheet->cid = 0;

    return sheet;
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

    if (sheet != NULL) {
	gdk_window_raise(sheet->win->window);
	return;
    }

    sheet = sheet_new();
    if (sheet == NULL) return;

    if (pdinfo == NULL && datainfo->v == 1) {
	errbox(_("Please add a variable to the dataset first"));
	return;
    }

    sheetwidth = get_obs_col_width() + 6 * get_data_col_width() + 14;

    sheet->win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(sheet->win), _("gretl: edit data"));
    gtk_window_set_default_size(GTK_WINDOW(sheet->win), sheetwidth, 400);

    main_vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 5); 
    gtk_container_add(GTK_CONTAINER(sheet->win), main_vbox);
    gtk_widget_show(main_vbox);

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
    gtk_widget_show(mbar);

    build_sheet_popup(&sheet->popup, sheet);

    status_box = gtk_hbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(status_box), 0);
    gtk_box_pack_start(GTK_BOX(main_vbox), status_box, FALSE, FALSE, 0);
    gtk_widget_show(status_box);

    sheet->locator = gtk_statusbar_new(); 
    gtk_widget_set_size_request(sheet->locator, get_obs_col_width(), 20);
    gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(sheet->locator), FALSE);
    sheet->cid = gtk_statusbar_get_context_id (GTK_STATUSBAR(sheet->locator), 
					       "current row and column");
    gtk_box_pack_start(GTK_BOX(status_box), sheet->locator, FALSE, FALSE, 0);
    gtk_widget_show(sheet->locator);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
                                    GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW(scroller), 
					 GTK_SHADOW_IN);

    gtk_box_pack_start(GTK_BOX(main_vbox), scroller, TRUE, TRUE, 0);
    gtk_widget_show(scroller);

    sheet->view = data_sheet_new(sheet, datainfo->n, datainfo->v);
    gtk_container_add(GTK_CONTAINER(scroller), sheet->view);
    gtk_widget_show(sheet->view);

    /* apply and close buttons */
    button_box = gtk_hbox_new (FALSE, 5);
    gtk_box_set_homogeneous (GTK_BOX (button_box), TRUE);
    gtk_box_pack_start (GTK_BOX (main_vbox), button_box, FALSE, FALSE, 0);
    gtk_widget_show(button_box);

    g_signal_connect(G_OBJECT(sheet->win), "destroy",
		     G_CALLBACK(free_spreadsheet), &sheet);

    tmp = gtk_button_new_from_stock(GTK_STOCK_APPLY);
    gtk_box_pack_start (GTK_BOX (button_box), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(get_data_from_sheet), sheet);
    gtk_widget_show(tmp);

    tmp = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_box_pack_start (GTK_BOX (button_box), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), sheet->win);
    gtk_widget_show(tmp);

    g_signal_connect (G_OBJECT(sheet->view), "button_press_event",
		      G_CALLBACK(popup_menu_handler),
		      sheet->popup);

    if (pdinfo != NULL)
	add_skel_to_sheet(sheet);
    else 
	add_data_to_sheet(sheet);

    gtk_widget_show(sheet->win);
}

