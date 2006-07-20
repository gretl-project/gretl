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
#include "dlgutils.h"
#include "menustate.h"

#include <errno.h>
#include <ctype.h>
#include <float.h>

#undef SSDEBUG
#undef CELLDEBUG

typedef enum {
    SHEET_SUBSAMPLED    = 1 << 1,
    SHEET_SHORT_VARLIST = 1 << 2,
    SHEET_INSERT_OBS_OK = 1 << 3,
    SHEET_ADD_OBS_OK    = 1 << 4
} SheetExtraFlags;

enum {
    SHEET_AT_END,
    SHEET_AT_POINT
};

typedef struct {
    GtkWidget *view;
    GtkWidget *win;
    GtkWidget *locator;
    GtkWidget *popup;
    GtkCellRenderer *dumbcell;
    GtkCellRenderer *datacell;
    gchar location[64];
    int *varlist;
    int datacols, datarows;
    int totcols;
    int orig_nobs;
    int added_vars;
    int orig_main_v;
    int modified;
    int flags;
    guint cid;
    guint point;
} Spreadsheet;

static void sheet_add_var_callback (gpointer data, guint code, GtkWidget *w);
static void sheet_add_obs_callback (gpointer data, guint where, GtkWidget *w);
static void get_data_from_sheet (GtkWidget *w, Spreadsheet *sheet);
static void set_up_sheet_column (GtkTreeViewColumn *column, gint width, 
				 gboolean expand);
static gint get_data_col_width (void);
static int add_data_column (Spreadsheet *sheet);
static void create_sheet_cell_renderers (Spreadsheet *sheet);
static void set_dataset_locked (gboolean s);

static GtkItemFactoryEntry sheet_items[] = {
    { N_("/_Observation"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Observation/_Append obs"), NULL, sheet_add_obs_callback, SHEET_AT_END, 
      NULL, GNULL },
    { N_("/Observation/_Insert obs"), NULL, sheet_add_obs_callback, SHEET_AT_POINT, 
      NULL, GNULL },
    { N_("/_Variable"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Variable/_Add"), NULL, sheet_add_var_callback, 0, NULL, GNULL }
};

static void disable_obs_menu (GtkItemFactory *ifac)
{
    GtkWidget *w = gtk_item_factory_get_item(ifac, "/Observation");

    gtk_widget_set_sensitive(w, FALSE);
}

static void disable_insert_obs_item (GtkItemFactory *ifac)
{
    GtkWidget *w = 
	gtk_item_factory_get_item(ifac, "/Observation/Insert obs");

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

static char *single_underscores (char *targ, const char *src)
{
    char *p = targ;

    while (*src) {
	if (*src == '_' && *(src + 1) == '_') {
	    src++;
	    *p++ = '_';
	} else {
	    *p++ = *src;
	}
	src++;
    }

    *p = '\0';

    return targ;
}

static void set_locator_label (Spreadsheet *sheet, GtkTreePath *path,
			       GtkTreeViewColumn *column)
{
    GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(sheet->view));
    GtkTreeIter iter;
    gchar *row_label;
    const gchar *col_label;
    char tmp[VNAMELEN];

    gtk_tree_model_get_iter(model, &iter, path);
    gtk_tree_model_get(model, &iter, 0, &row_label, -1);
    col_label = gtk_tree_view_column_get_title(column);
    single_underscores(tmp, col_label);
    sprintf(sheet->location, "%s, %s", tmp, row_label);
    gtk_statusbar_pop(GTK_STATUSBAR(sheet->locator), sheet->cid);
    gtk_statusbar_push(GTK_STATUSBAR(sheet->locator), 
		       sheet->cid, sheet->location);
    g_free(row_label);
}

static void move_to_next_cell (Spreadsheet *sheet, GtkTreePath *path,
			       GtkTreeViewColumn *column)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    gint nextrow, colnum;

    nextrow = gtk_tree_path_get_indices(path)[0] + 1;

    /* go to next row, if possible... */
    if (nextrow < sheet->datarows) {
	GtkTreePath *newpath;
	gchar pstr[8];

	sprintf(pstr, "%d", nextrow);
	newpath = gtk_tree_path_new_from_string(pstr);
	if (newpath != NULL) {
	    gtk_tree_view_set_cursor(view, newpath, column, FALSE);
	    set_locator_label(sheet, newpath, column);
	    gtk_tree_path_free(newpath);
	}
    } else {
	/* ...or try the next column */
	gpointer p = g_object_get_data(G_OBJECT(column), "colnum");

	if (p != NULL && (colnum = GPOINTER_TO_INT(p)) < sheet->datacols) {
	    GtkTreeViewColumn *nextcol = gtk_tree_view_get_column(view, colnum + 1);

	    if (nextcol != NULL) {
		gtk_tree_view_set_cursor(view, path, nextcol, FALSE);
		set_locator_label(sheet, path, nextcol);
	    }		
	}
    }
    /* couldn't find a "next cell" to go to */
}

static gint sheet_cell_edited (GtkCellRendererText *cell,
			       const gchar *path_string,
			       const gchar *user_text,
			       Spreadsheet *sheet)
{
    const gchar *new_text = NULL;
    int err = 0;

    if (!strcmp(user_text, "na") || !strcmp(user_text, "NA")) {
	/* allow conversion to "missing" */
	new_text = "";
    } else {
	err = check_atof(user_text);
	if (err) {
	    errbox(get_gretl_errmsg());
	} else {
	    new_text = user_text;
	}
    }

    if (!err && new_text != NULL) {
	GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
	GtkTreeModel *model = gtk_tree_view_get_model(view);
	GtkTreeViewColumn *column;
	GtkTreePath *path;
	GtkTreeIter iter;
	gchar *old_text;
	gint colnum;

	gtk_tree_view_get_cursor(view, NULL, &column);
	colnum = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(column), "colnum"));
	path = gtk_tree_path_new_from_string(path_string);
	gtk_tree_model_get_iter(model, &iter, path);
	gtk_tree_model_get(model, &iter, colnum, &old_text, -1);

	if (old_text != NULL && strcmp(old_text, new_text)) {
	    gtk_list_store_set(GTK_LIST_STORE(model), &iter, 
			       colnum, new_text, -1);
	    sheet->modified = 1;
	}
	move_to_next_cell(sheet, path, column);
	gtk_tree_path_free(path);
	g_free(old_text);
    }

    return FALSE;
}

static void 
spreadsheet_scroll_to_new_col (Spreadsheet *sheet, GtkTreeViewColumn *column)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreePath *path;
    gchar *pathstr;
    GtkAdjustment *adj;
    GtkWidget *vp;

    pathstr = g_strdup("0");
    path = gtk_tree_path_new_from_string(pathstr);
    gtk_tree_view_set_cursor(view, path, column, TRUE);
    vp = gtk_widget_get_ancestor(GTK_WIDGET(view), GTK_TYPE_BIN);
    adj = gtk_viewport_get_hadjustment(GTK_VIEWPORT(vp));
    gtk_adjustment_set_value(adj, adj->upper);
    gtk_tree_path_free(path);
    g_free(pathstr);
}

static int real_add_new_var (Spreadsheet *sheet, const char *varname)
{
    GtkTreeViewColumn *column;
    gint cols, colnum;
    char tmp[32];

    if (add_data_column(sheet)) {
	return 1;
    }

#if SSDEBUG
    fprintf(stderr, "real_add_new_var, after add_data_column: sheet->totcols=%d\n", 
	    sheet->totcols);
#endif

    column = gtk_tree_view_column_new();
    double_underscores(tmp, varname);
    gtk_tree_view_column_set_title(column, tmp);
    set_up_sheet_column(column, get_data_col_width(), TRUE);

    cols = gtk_tree_view_insert_column(GTK_TREE_VIEW(sheet->view), column, -1);
    colnum = cols - 1;

#if 0	
    create_sheet_cell_renderers(sheet);
#endif

    gtk_tree_view_column_pack_start(column, sheet->datacell, TRUE);
    gtk_tree_view_column_set_attributes(column, 
					sheet->datacell,
					"text", colnum,
					NULL);
    g_object_set_data(G_OBJECT(column), "colnum", GINT_TO_POINTER(colnum));

    /* scroll to editing position if need be */
    spreadsheet_scroll_to_new_col(sheet, column);

    sheet->modified = 1;

    return 0;
}

static void 
spreadsheet_scroll_to_foot (Spreadsheet *sheet, int row, int col)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreePath *path;
    gchar *pathstr;
    GtkTreeViewColumn *column;
    GtkAdjustment *adj;
    GtkWidget *vp;

    pathstr = g_strdup_printf("%d", row);
    path = gtk_tree_path_new_from_string(pathstr);
    column = gtk_tree_view_get_column(view, col);
    gtk_tree_view_set_cursor(view, path, column, FALSE);
    vp = gtk_widget_get_ancestor(GTK_WIDGET(view), GTK_TYPE_BIN);
    adj = gtk_viewport_get_vadjustment(GTK_VIEWPORT(vp));
    gtk_adjustment_set_value(adj, adj->upper);
    gtk_tree_path_free(path);
    g_free(pathstr);
}

static void 
real_add_new_obs (Spreadsheet *sheet, const char *obsname, int n)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    gint rownum = 0;
    gint oldrows = sheet->datarows;
    GtkListStore *store;
    gchar *pathstr = NULL;
    GtkTreeIter iter;
    gchar rowlabel[10];
    gint i, j;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));

    if (sheet->point == SHEET_AT_END) {
	rownum = sheet->datarows - 1;
	pathstr = g_strdup_printf("%d", rownum);
	for (i=0; i<n; i++) {
	    gtk_list_store_append(store, &iter);
	}
    } else if (sheet->point == SHEET_AT_POINT) {
	GtkTreePath *path;

	gtk_tree_view_get_cursor(view, &path, NULL);
	rownum = gtk_tree_path_get_indices(path)[0];
	gtk_tree_model_get_iter(GTK_TREE_MODEL(store), &iter, path);
	gtk_list_store_insert(store, &iter, rownum);
	gtk_tree_path_free(path);
    } else {
	return;
    }

    if (datainfo->markers && obsname != NULL) {
	gtk_list_store_set(store, &iter, 0, obsname, -1);
    } else if (sheet->point == SHEET_AT_END) {
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter,
					    pathstr);
	for (j=0; j<n; j++) {
	    gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	    get_full_obs_string(rowlabel, sheet->datarows + j, datainfo);
	    gtk_list_store_set(store, &iter, 0, rowlabel, -1);
	}
    }

    sheet->datarows += n;

    if (pathstr != NULL) {
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter,
					    pathstr);
	for (j=0; j<n; j++) {
	    gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	    for (i=1; i<=sheet->datacols; i++) {
		gtk_list_store_set(store, &iter, i, "", -1);
	    }	
	}
    } else {
	for (i=1; i<=sheet->datacols; i++) {
	    gtk_list_store_set(store, &iter, i, "", -1);
	}
    }	

    if (sheet->point == SHEET_AT_POINT && !datainfo->markers) {
	for (i=rownum; i<sheet->datarows; i++) {
	    get_full_obs_string(rowlabel, i, datainfo);
	    gtk_list_store_set(store, &iter, 0, rowlabel, -1);
	    gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	}
    } 

    if (sheet->point == SHEET_AT_END) {
	spreadsheet_scroll_to_foot(sheet, oldrows, 1);
    } else {
	GtkTreePath *path;
	GtkTreeIter insiter;

	pathstr = g_strdup_printf("%d", rownum);
	path = gtk_tree_path_new_from_string(pathstr);
	gtk_tree_model_get_iter(GTK_TREE_MODEL(store), &insiter, path);
	gtk_list_store_set(store, &insiter, 1, "", -1);
	gtk_tree_path_free(path);
    }

    if (pathstr != NULL) {
	g_free(pathstr);
    }

    sheet->modified = 1;
}

static void name_new_var (GtkWidget *widget, dialog_t *dlg) 
{
    Spreadsheet *sheet = (Spreadsheet *) edit_dialog_get_data(dlg);
    const gchar *buf;
    char varname[VNAMELEN];

    buf = edit_dialog_get_text(dlg);

    if (buf == NULL || validate_varname(buf)) return;

    *varname = 0;
    strncat(varname, buf, VNAMELEN - 1);

    close_dialog(dlg);

    if (real_add_new_var(sheet, varname)) {
	errbox(_("Out of memory attempting to add variable"));
    }
}

static void name_new_obs (GtkWidget *widget, dialog_t *dlg) 
{
    Spreadsheet *sheet = (Spreadsheet *) edit_dialog_get_data(dlg);
    const gchar *buf;
    char obsmarker[OBSLEN];

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    *obsmarker = 0;
    strncat(obsmarker, buf, OBSLEN - 1);

    close_dialog(dlg);
    real_add_new_obs(sheet, obsmarker, 1);
}

static void name_var_dialog (Spreadsheet *sheet) 
{
    edit_dialog (_("gretl: name variable"), 
		 _("Enter name for new variable\n"
		   "(max. 15 characters)"),
		 NULL, name_new_var, sheet, 
		 0, VARCLICK_NONE, NULL);
}

static void new_case_dialog (Spreadsheet *sheet) 
{
    edit_dialog (_("gretl: case marker"), 
		 _("Enter case marker for new obs\n"
		   "(max. 8 characters)"),
		 NULL, name_new_obs, sheet, 
		 0, VARCLICK_NONE, NULL);
}

static int add_data_column (Spreadsheet *sheet)
{
    GType *types;
    GtkListStore *old_store, *new_store;
    GtkTreeIter old_iter, new_iter;
    gint i, row;

    /* This is relatively complex because, so far as I can tell, you can't
       append or insert additional columns in a GtkListStore: we have to
       create a whole new liststore and copy the old info across.
    */

    /* make an expanded column types list */
    types = mymalloc((sheet->totcols + 1) * sizeof *types);
    if (types == NULL) {
	return 1;
    }

    sheet->datacols += 1;
    sheet->totcols += 1;

    /* configure the types */
    for (i=0; i<sheet->totcols; i++) {
	types[i] = G_TYPE_STRING;
    }

    /* get pointers to original and new stores */
    old_store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sheet->view)));
    new_store = gtk_list_store_newv(sheet->totcols, types);
    free(types);

    /* go to start of old and new lists */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(old_store), &old_iter);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(new_store), &new_iter);

    /* construct the new table */
    for (row=0; row<sheet->datarows; row++) {
	gchar *str;
	int col;

	gtk_list_store_append(new_store, &new_iter);

	for (col=0; col<sheet->totcols; col++) {
	    if (col < sheet->datacols) {
		/* copy labels and original data */
		gtk_tree_model_get(GTK_TREE_MODEL(old_store), &old_iter, col, &str, -1);
		gtk_list_store_set(new_store, &new_iter, col, str, -1);
		g_free(str);
	    } else { 
		/* new data values: set blank */
		gtk_list_store_set(new_store, &new_iter, col, "", -1);
	    }
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(old_store), &old_iter);
    }    

    gtk_tree_view_set_model(GTK_TREE_VIEW(sheet->view), GTK_TREE_MODEL(new_store));
    g_object_unref(G_OBJECT(new_store));

    sheet->added_vars += 1;

#if SSDEBUG
    fprintf(stderr, "add_data_column: sheet->added_vars now = %d\n",
	    sheet->added_vars);
#endif

    return 0;
}

static void sheet_add_var_callback (gpointer data, guint u, GtkWidget *w)
{
    Spreadsheet *sheet = (Spreadsheet *) data;

    name_var_dialog(sheet);
}

static void sheet_add_obs_callback (gpointer data, guint where, GtkWidget *w)
{
    Spreadsheet *sheet = (Spreadsheet *) data;

    sheet->point = where;

    if (datainfo->markers) {
	new_case_dialog(sheet);
    } else if (where == SHEET_AT_END) {
	int n = add_obs_dialog(NULL, 1);

	if (n > 0) {
	    real_add_new_obs(sheet, NULL, n);
	}
    } else {
	real_add_new_obs(sheet, NULL, 1);
    }
}

static void popup_sheet_add_obs (GtkWidget *w, Spreadsheet *sheet)
{
    sheet_add_obs_callback(sheet, SHEET_AT_END, NULL);
}

static void popup_sheet_insert_obs (GtkWidget *w, Spreadsheet *sheet)
{
    sheet_add_obs_callback(sheet, SHEET_AT_POINT, NULL);
}

static void popup_sheet_add_var (GtkWidget *w, Spreadsheet *sheet)
{
    sheet_add_var_callback(sheet, 0, NULL);
}

static void build_sheet_popup (Spreadsheet *sheet)
{
    if (sheet->popup != NULL) return;

    sheet->popup = gtk_menu_new();

    add_popup_item(_("Add Variable"), sheet->popup, 
		   G_CALLBACK(popup_sheet_add_var),
		   sheet);

    if (sheet->flags & SHEET_ADD_OBS_OK) {
	add_popup_item(_("Add Observation"), sheet->popup,
		       G_CALLBACK(popup_sheet_add_obs),
		       sheet);
    } 

    if (sheet->flags & SHEET_INSERT_OBS_OK) {
	add_popup_item(_("Insert Observation"), sheet->popup,
		       G_CALLBACK(popup_sheet_insert_obs),
		       sheet);
    }
}

static gboolean update_cell_position (GtkTreeView *view, Spreadsheet *sheet)
{
    GtkTreePath *path = NULL;
    GtkTreeViewColumn *column;
    static gint oldrow, oldcol;

    /* this is connected to the "cursor-changed" signal */

#if CELLDEBUG
    fprintf(stderr, "** update_cell_position()\n");
#endif

    gtk_tree_view_get_cursor(view, &path, &column);

    if (path != NULL && column != NULL) {
	gint newrow = gtk_tree_path_get_indices(path)[0];
	gint newcol = 
	    GPOINTER_TO_INT(g_object_get_data(G_OBJECT(column), "colnum"));

	if (newcol == 0) {
	    /* not a data column */
	    gtk_tree_path_free(path);
	    return TRUE;
	}

	if (newrow != oldrow || newcol != oldcol) {
#if CELLDEBUG
	    fprintf(stderr, " activating cell(%d, %d)\n", newrow, newcol);
#endif
	    set_locator_label(sheet, path, column);
	    oldrow = newrow;
	    oldcol = newcol;
	    gtk_tree_view_set_cursor(view, path, column, 
				     FALSE); /* start editing? */
	} else {
#if CELLDEBUG
	   fprintf(stderr, " still in cell(%d, %d)\n", oldrow, oldcol); 
#endif
	}
    }

    if (path != NULL) {
	gtk_tree_path_free(path);
    }

    return TRUE; /* is this right? */
}

static int 
var_added_since_ssheet_opened (int i, Spreadsheet *sheet, int main_v)
{
    return (i >= sheet->orig_main_v && i < main_v);
}

/* pull modified values from the data-editing spreadsheet
   into the main dataset */

static void get_data_from_sheet (GtkWidget *w, Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkTreeViewColumn *column;
    GtkTreeModel *model;
    int oldv = datainfo->v;
    int newvars = sheet->added_vars;
    int newobs = sheet->datarows - sheet->orig_nobs;
    int missobs = 0;
    int i, colnum, s, t;

    if (!sheet->modified) {
	infobox(_("No changes were made"));
	return;
    }    

#if SSDEBUG
    fprintf(stderr, "get_data_from_sheet: oldv = %d, newvars = %d, newobs = %d\n",
	    oldv, newvars, newobs);
#endif

    model = gtk_tree_view_get_model(view);

    /* first extend series length, if needed */

    if (newobs > 0) {
	if (dataset_add_observations(newobs, &Z, datainfo, OPT_A) ||
	    dataset_destroy_hidden_variables(&Z, datainfo, 0)) {
	    nomem();
	    return;
	}
    }

    /* then add any new variables to data set */

    if (newvars > 0) {
	int vi, v0 = sheet->varlist[0];
	const gchar *newname;

	if (dataset_add_series(newvars, &Z, datainfo)) {
	    nomem();
	    return;
	}

	for (i=0; i<newvars; i++) { 
	    vi = oldv + i;
	    gretl_list_append_term(&sheet->varlist, vi);
	    if (sheet->varlist == NULL) {
		nomem();
		return;
	    }
	    colnum = v0 + 1 + i;
	    column = gtk_tree_view_get_column(view, colnum);
	    newname = gtk_tree_view_column_get_title(column);
	    strcpy(datainfo->varname[vi], newname); 
#if SSDEBUG
	    fprintf(stderr, " added var %d (%s) from column %d\n",
		    vi, newname, colnum);
#endif
	}
    }

    /* copy data values from spreadsheet */

    colnum = 1;
    for (i=1; i<=sheet->varlist[0]; i++) {
	int vi = sheet->varlist[i];
	gchar *numstr;

	gtk_tree_model_get_iter_first(model, &iter);

#if SSDEBUG
	fprintf(stderr, " updating data for var %d (%s) from column %d\n",
		vi, datainfo->varname[vi], colnum);
#endif
	t = datainfo->t1;
	for (s=0; s<sheet->datarows; s++) {
	    gtk_tree_model_get(model, &iter, colnum, &numstr, -1);
	    if (*numstr != '\0') {
		Z[vi][t++] = atof(numstr); 
	    } else {
		Z[vi][t++] = NADBL;
		missobs = 1;
	    }
	    g_free(numstr);
	    gtk_tree_model_iter_next(model, &iter);
	}
	colnum++;
    }

    /* copy observation markers, if relevant */

    if (datainfo->markers && datainfo->S != NULL) {
	gchar *marker;

	gtk_tree_model_get_iter_first(model, &iter);
	t = datainfo->t1;
	for (s=0; s<sheet->datarows; s++) {
	    gtk_tree_model_get(model, &iter, 0, &marker, -1);
	    strcpy(datainfo->S[t++], marker);
	    g_free(marker);
	    gtk_tree_model_iter_next(model, &iter);
	}
    }

    data_status |= (GUI_DATA | MODIFIED_DATA);
    register_data(NULL, NULL, 0);

    if (missobs) {
	infobox(_("Warning: there were missing observations"));
    } 

    sheet->modified = 0;
}

static void select_first_editable_cell (Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreePath *path;
    GtkTreeViewColumn *column;

    path = gtk_tree_path_new_from_string("0");
    column = gtk_tree_view_get_column(view, 1);
    gtk_tree_view_set_cursor(view, path, column, FALSE);
    set_locator_label(sheet, path, column);

    gtk_tree_path_free(path);
}

static int add_data_to_sheet (Spreadsheet *sheet, SheetCode code)
{
    gchar rowlabel[OBSLEN];
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkListStore *store;
    int i, t;

#if SSDEBUG
    fprintf(stderr, "Doing add_data_to_sheet\n");
#endif

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));

    /* insert observation markers */

    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	get_full_obs_string(rowlabel, t, datainfo);
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, rowlabel, -1);
    }

    sheet->datarows = datainfo->t2 - datainfo->t1 + 1;

    /* now insert data values */

    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (code == SHEET_NEW_DATASET) {
	    /* no hidden vars to consider; insert NAs for first var */
	    gtk_list_store_set(store, &iter, 1, "", -1);
	} else {
	    char numstr[32];
	    int colnum = 0;
	    int vi;

	    for (i=1; i<=sheet->varlist[0]; i++) {
		vi = sheet->varlist[i];
		if (na(Z[vi][t])) {
		    *numstr = '\0';
		} else {
		    sprintf(numstr, "%.*g", DBL_DIG, Z[vi][t]);
		}
		gtk_list_store_set(store, &iter, ++colnum, numstr, -1);
	    }
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }

    sheet->orig_main_v = datainfo->v;

#if SSDEBUG
    fprintf(stderr, " datarows=%d, orig_vars=%d, orig_main_v=%d\n", 
	    sheet->datarows, sheet->varlist[0], sheet->orig_main_v);
#endif

    return 0;
}

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

    if (width == 0) {
	width = get_string_width("XXXXXXXXX");
    }
    return width;
}

static gint get_data_col_width (void)
{
    static gint width;

    if (width == 0) {
	width = get_string_width("-00.000000");
    }
    return width;
}

static void 
set_up_sheet_column (GtkTreeViewColumn *column, gint width, gboolean expand)
{
    gtk_tree_view_column_set_alignment(column, 0.5); /* header centered */
    gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_AUTOSIZE);
    gtk_tree_view_column_set_resizable(column, TRUE);
    gtk_tree_view_column_set_min_width(column, width);
    gtk_tree_view_column_set_expand(column, expand);
}

static void create_sheet_cell_renderers (Spreadsheet *sheet)
{
    GtkCellRenderer *r;

    r = gtk_cell_renderer_text_new();
    g_object_set(r, "ypad", 1, 
		 "xalign", 1.0,
		 "background", "gray",
		 "editable", FALSE, NULL);
    sheet->dumbcell = r;

    r = gtk_cell_renderer_text_new();
    g_object_set(r, "ypad", 1, 
		 "xalign", 1.0, 
		 "editable", TRUE, NULL);
    g_signal_connect (G_OBJECT (r), "edited",
		      G_CALLBACK(sheet_cell_edited), sheet);
    sheet->datacell = r;
}

static void manufacture_keystroke (GtkWidget *widget, guint uval)
{
    GdkKeymapKey *keys;
    gint n_keys;

    if (gdk_keymap_get_entries_for_keyval(NULL, uval, &keys, &n_keys)) {
	guint16 hardware_keycode;
	GdkEvent *event;

	hardware_keycode = keys[0].keycode;
	g_free(keys);

	event = gdk_event_new(GDK_KEY_PRESS);
	event->key.window = g_object_ref(widget->window);
	event->key.hardware_keycode = hardware_keycode;

	event->key.keyval = gdk_unicode_to_keyval(uval);
	event->key.length = 1;

	event->key.send_event = FALSE;
	event->key.time = GDK_CURRENT_TIME;   

	gtk_main_do_event(event);
	gdk_event_free(event);
    }
}

static gint catch_spreadsheet_key (GtkWidget *view, GdkEventKey *key, 
				   Spreadsheet *sheet)
{
    if (key->keyval == GDK_Tab) {
	/* FIXME */
	;
    }

    /* prevent cursor movement outside of the data area */

    if (key->keyval == GDK_Right || key->keyval == GDK_Left) {
	GtkTreeViewColumn *column;
	gpointer p;

	gtk_tree_view_get_cursor(GTK_TREE_VIEW(view), NULL, &column);
	p = g_object_get_data(G_OBJECT(column), "colnum");
	if (p != NULL) {
	    int colnum = GPOINTER_TO_INT(p);

	    if (key->keyval == GDK_Left && colnum == 1) {
		return TRUE;
	    }
	}
    } 

    /* numeric key: start editing */

    else if ((key->keyval >= GDK_0 && key->keyval <= GDK_9) ||
	key->keyval == GDK_minus || key->keyval == GDK_period) {
	GtkTreePath *path = NULL;
	GtkTreeViewColumn *column;

	gtk_tree_view_get_cursor(GTK_TREE_VIEW(view), &path, &column);

	if (path != NULL && column != NULL) {
	    gtk_tree_view_set_cursor(GTK_TREE_VIEW(view), path, column, 
				     TRUE);
	    gtk_tree_path_free(path);
	    manufacture_keystroke(view, key->keyval);
	}
    }

    return FALSE;
}

/* fuller version of "catch_spreadsheet_click" for gtk 2.2 and higher */

static gint catch_spreadsheet_click (GtkWidget *view, GdkEvent *event,
				     Spreadsheet *sheet)
{   
    GdkModifierType mods; 
    gint ret = FALSE;

# if CELLDEBUG
    fprintf(stderr, "** catch_spreadsheet_click()\n");
# endif

    if (event->type != GDK_BUTTON_PRESS) {
	return FALSE;
    }

    gdk_window_get_pointer(view->window, NULL, NULL, &mods);

    if (mods & GDK_BUTTON3_MASK) {
	GdkEventButton *bevent = (GdkEventButton *) event;

	if (sheet->popup == NULL)
	    build_sheet_popup(sheet);

	gtk_menu_popup (GTK_MENU(sheet->popup), NULL, NULL, NULL, NULL,
			bevent->button, bevent->time);
	return TRUE;
    }	    
	
    if (mods & GDK_BUTTON1_MASK) {
	GdkEventButton *bevent = (GdkEventButton *) event;
	GtkTreePath *path;
	GtkTreeViewColumn *column;

# if CELLDEBUG
	fprintf(stderr, "Got button 1 click\n");
# endif

	gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(sheet->view),
				      (gint) bevent->x, 
				      (gint) bevent->y,
				      &path, &column,
				      NULL, NULL);
	if (column != NULL) {
	    gpointer p = g_object_get_data(G_OBJECT(column), "colnum");
	    gint colnum = GPOINTER_TO_INT(p);

# if CELLDEBUG
	    fprintf(stderr, "Clicked column: colnum = %d\n", colnum);
# endif

	    if (colnum == 0) {
		/* don't respond to a click in a non-data column */
		ret = TRUE;
	    } else {
		/* otherwise start editing on clicked cell */
		gtk_tree_view_set_cursor(GTK_TREE_VIEW(sheet->view), 
					 path, column, TRUE);
		ret = TRUE;
	    }
	}
	gtk_tree_path_free(path);
    }

# if CELLDEBUG
    fprintf(stderr, "catch_spreadsheet_click returning %d\n", ret);
# endif

    return ret;
}

static int build_sheet_view (Spreadsheet *sheet)
{
    GtkListStore *store; 
    GtkWidget *view;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    GType *types;
    gint i, width, colnum;

    sheet->datacols = 0;

    for (i=1; i<=sheet->varlist[0]; i++) {
	if (spreadsheet_hide(sheet->varlist[i], datainfo)) {
	    gretl_list_delete_at_pos(sheet->varlist, i--);
	} else {
	    sheet->datacols += 1;
	}
    }  

    if (sheet->varlist[0] < datainfo->v - 1) {
	sheet->flags |= SHEET_SHORT_VARLIST;
    }

    /* obs col, data cols */
    sheet->totcols = sheet->datacols + 1;

    types = mymalloc(sheet->totcols * sizeof *types);
    if (types == NULL) {
	return 1;
    }

    /* configure the types */
    for (i=0; i<sheet->totcols; i++) {
	types[i] = G_TYPE_STRING;
    }

    store = gtk_list_store_newv(sheet->totcols, types);
    free(types);

    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), TRUE);

    /* build and attach the (two) cell renderers */
    create_sheet_cell_renderers(sheet);

    /* construct the observation marker column */
    width = get_obs_col_width();
    column = gtk_tree_view_column_new_with_attributes(NULL,
						      sheet->dumbcell,
						      "text", 0, 
						      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
    set_up_sheet_column(column, width, FALSE);
    g_object_set_data(G_OBJECT(column), "colnum", GINT_TO_POINTER(0));

    /* construct the data columns */
    width = get_data_col_width();
    colnum = 0;
    for (i=1; i<=sheet->varlist[0]; i++) {
	int vi = sheet->varlist[i];
	char tmp[16];

	if (spreadsheet_hide(vi, datainfo)) {
	    continue;
	}
	colnum++;
	double_underscores(tmp, datainfo->varname[vi]);
	column = gtk_tree_view_column_new_with_attributes(tmp,
							  sheet->datacell,
							  "text", 
							  colnum, 
							  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	set_up_sheet_column(column, width, TRUE);
	g_object_set_data(G_OBJECT(column), "colnum", GINT_TO_POINTER(colnum));
    }

    /* set the selection property on the tree view */
    select = gtk_tree_view_get_selection (GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_NONE);

    g_signal_connect (G_OBJECT(view), "cursor-changed",
		      G_CALLBACK(update_cell_position), sheet);
    g_signal_connect (G_OBJECT(view), "key_press_event",
		      G_CALLBACK(catch_spreadsheet_key), sheet);

    /* attach to sheet struct */
    sheet->view = view;

    return 0;
}

static void free_spreadsheet (GtkWidget *widget, Spreadsheet **psheet) 
{
    Spreadsheet *sheet = *psheet;

    if (sheet->popup != NULL) {
	gtk_widget_destroy(sheet->popup);
    }

    free(sheet->varlist);

    free(sheet);
    *psheet = NULL;

    set_dataset_locked(FALSE);
}

static Spreadsheet *spreadsheet_new (int flags)
{
    Spreadsheet *sheet;

    sheet = mymalloc(sizeof *sheet);
    if (sheet == NULL) return NULL;

    sheet->view = NULL;
    sheet->win = NULL;
    sheet->locator = NULL;
    sheet->popup = NULL;
    sheet->dumbcell = NULL;
    sheet->datacell = NULL;
    sheet->datacols = sheet->datarows = 0;
    sheet->totcols = 0;
    sheet->added_vars = 0;
    sheet->orig_main_v = 0;
    sheet->modified = 0;
    sheet->cid = 0;
    sheet->flags = flags;

    if (datainfo->t1 != 0 || datainfo->t2 < datainfo->n - 1) {
	sheet->flags |= SHEET_SUBSAMPLED;
    }

    sheet->orig_nobs = datainfo->t2 - datainfo->t1 + 1;

    if (sheet->flags & SHEET_NEW_DATASET) {
	sheet->varlist = gretl_list_new(1);
    } else {
	sheet->varlist = main_window_selection_as_list();
    }    

    if (sheet->varlist == NULL) {
	free(sheet);
	sheet = NULL;
    } else if (sheet->flags & SHEET_NEW_DATASET) {
	sheet->varlist[1] = 1;
    }

    return sheet;
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

static gint maybe_exit_sheet (GtkWidget *w, Spreadsheet *sheet)
{
    int resp;

    if (sheet->modified) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    get_data_from_sheet(NULL, sheet);
	} else if (resp == GRETL_CANCEL || resp == -1) {
	    return FALSE;
	}
    } 

    if (sheet->flags & SHEET_NEW_DATASET) {
	empty_dataset_guard();
    }
  
    gtk_widget_destroy(sheet->win);

    return FALSE;
}

static gint sheet_delete_event (GtkWidget *w, GdkEvent *event,
				Spreadsheet *sheet)
{
    int resp;

    if (sheet->modified) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    get_data_from_sheet(NULL, sheet);
	} else if (resp == GRETL_CANCEL) {
	    return TRUE;
	}
    }

    if (sheet->flags & SHEET_NEW_DATASET) {
	empty_dataset_guard();
    }

    return FALSE;
}

static void 
sheet_adjust_menu_state (Spreadsheet *sheet, GtkItemFactory *ifac)
{
    sheet->flags |= SHEET_ADD_OBS_OK | SHEET_INSERT_OBS_OK;

    if (complex_subsampled() || datainfo->t2 < datainfo->n - 1) {
	sheet->flags &= ~SHEET_ADD_OBS_OK;
	sheet->flags &= ~SHEET_INSERT_OBS_OK;
	disable_obs_menu(ifac);
    } else if (sheet->flags & SHEET_SHORT_VARLIST) {
	sheet->flags &= ~SHEET_INSERT_OBS_OK;
	disable_insert_obs_item(ifac);
    }
}

void show_spreadsheet (SheetCode code) 
{
    static Spreadsheet *sheet;    
    GtkWidget *tmp, *button_box;
    GtkWidget *scroller, *main_vbox;
    GtkWidget *hbox, *padbox;
    GtkWidget *status_box, *mbar;
    GtkItemFactory *ifac;
    int sheetwidth;
    int err = 0;

#ifdef G_OS_WIN32
    if (datainfo->t2 - datainfo->t1 > 1600) {
	errbox(_("Sorry, can't edit more than 1600 rows"));
	return;
    }    
#endif

    if (sheet != NULL) {
	gdk_window_raise(sheet->win->window);
	return;
    }

    sheet = spreadsheet_new(code);
    if (sheet == NULL) {
	return;
    }

    if (datainfo->v == 1) {
	errbox(_("Please add a variable to the dataset first"));
	return;
    }

    sheetwidth = get_obs_col_width() + 6 * get_data_col_width() + 14;

    sheet->win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(sheet->win), _("gretl: edit data"));
    gtk_window_set_default_size(GTK_WINDOW(sheet->win), sheetwidth, 400);

    g_signal_connect(G_OBJECT(sheet->win), "delete_event",
		     G_CALLBACK(sheet_delete_event), sheet);

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

    status_box = gtk_hbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(status_box), 0);
    gtk_box_pack_start(GTK_BOX(main_vbox), status_box, FALSE, FALSE, 0);
    gtk_widget_show(status_box);

    sheet->locator = gtk_statusbar_new(); 
    gtk_widget_set_size_request(sheet->locator, 2 * get_obs_col_width(), 20);
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

    hbox = gtk_hbox_new(FALSE, 1);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), hbox);
    gtk_widget_show(hbox);

    build_sheet_view(sheet);
    gtk_container_add(GTK_CONTAINER(hbox), sheet->view);
    gtk_widget_show(sheet->view); 

    sheet_adjust_menu_state(sheet, ifac);

    if (sheet->varlist[0] < 7) {
	/* padding to fill out the spreadsheet */
	padbox = gtk_vbox_new(FALSE, 1);
	gtk_container_add(GTK_CONTAINER(hbox), padbox);
	gtk_widget_show(padbox);
    }

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
    gtk_box_pack_start(GTK_BOX (button_box), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(maybe_exit_sheet), sheet);
    gtk_widget_show(tmp);

    g_signal_connect(G_OBJECT(sheet->view), "button_press_event",
		     G_CALLBACK(catch_spreadsheet_click),
		     sheet);

    add_data_to_sheet(sheet, code);

    if (err) {
	gtk_widget_destroy(sheet->win);
	return;
    }

    select_first_editable_cell(sheet);

    gtk_widget_show(sheet->win);

    /* we can't have the user making confounding changes elsewhere,
       while editing the dataset here */
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


