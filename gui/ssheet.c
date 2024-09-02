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

#include "gretl.h"
#include "gui_utils.h"
#include "treeutils.h"
#include "ssheet.h"
#include "dlgutils.h"
#include "menustate.h"
#include "selector.h"
#include "session.h"
#include "winstack.h"
#include "toolbar.h"
#include "cmdstack.h"
#include "fncall.h"
#include "usermat.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "gretl_cmatrix.h"

#include <errno.h>
#include <ctype.h>
#include <float.h>

#define SSDEBUG 0
#define CELLDEBUG 0

typedef enum {
    SHEET_SUBSAMPLED    = 1 << 0,
    SHEET_SHORT_VARLIST = 1 << 1,
    SHEET_INSERT_OBS_OK = 1 << 2,
    SHEET_ADD_OBS_OK    = 1 << 3,
    SHEET_USE_COMMA     = 1 << 4,
    SHEET_MODIFIED      = 1 << 5,
    SHEET_COLNAME_MOD   = 1 << 6,
    SHEET_CUSTOM_FMT    = 1 << 7
} SheetFlags;

enum {
    SHEET_AT_END,
    SHEET_AT_POINT
};

enum {
    NEXT_DOWN,
    NEXT_RIGHT,
    NEXT_UP,
    NEXT_LEFT,
    NEXT_SAME
};

typedef struct {
    GtkWidget *view;
    GtkWidget *win;
    GtkWidget *locator;
    GtkWidget *entry;
    GtkWidget *popup;
    GtkWidget *save;
    GtkWidget *apply;
    GtkUIManager *ui;
    GtkCellRenderer *textcell;
    GtkCellRenderer *datacell;
    GdkPixbuf *pbuf;
    gchar location[64];
    gretl_matrix *matrix;
    gretl_matrix *oldmat;
    char mname[VNAMELEN];
    const char **colnames;
    const char **rownames;
    int *varlist;
    int datacols, datarows;
    int totcols;
    int orig_nobs;
    int added_vars;
    int orig_main_v;
    int next;
    SheetCmd cmd;
    SheetFlags flags;
    guint cid;
    guint point;
    int *edits;
    int *inserts;
    int digits;
    char numfmt[8];
} Spreadsheet;

#define editing_series(s) (s->cmd == SHEET_EDIT_VARLIST || \
                           s->cmd == SHEET_EDIT_DATASET || \
                           s->cmd == SHEET_NEW_DATASET)

#define SERIES_DIGITS 10

static void set_up_sheet_column (GtkTreeViewColumn *column, gint width,
				 gboolean expand);
static gint get_data_col_width (void);
static int add_data_column (Spreadsheet *sheet);
static void create_sheet_cell_renderers (Spreadsheet *sheet);

static void matrix_fill_callback (GtkAction *action, gpointer data);
static void matrix_props_callback (GtkAction *action, gpointer data);
static void matrix_format_callback (GtkAction *action, gpointer data);
static void matrix_edit_callback (GtkAction *action, gpointer data);
static int update_sheet_from_matrix (Spreadsheet *sheet);
static void size_matrix_window (Spreadsheet *sheet);
static void set_ok_transforms (Spreadsheet *sheet);

static void sheet_show_popup (GtkWidget *w, Spreadsheet *sheet);
static void get_data_from_sheet (GtkWidget *w, Spreadsheet *sheet);
static gint maybe_exit_sheet (GtkWidget *w, Spreadsheet *sheet);

static void add_scalar_callback (GtkWidget *w, Spreadsheet *sheet);
static void update_scalars_from_sheet (Spreadsheet *sheet);
static void scalars_changed_callback (void);

static void
spreadsheet_scroll_to_foot (Spreadsheet *sheet, int row, int col);

enum {
    SHEET_ADD_BTN,
    SHEET_APPLY_BTN
};

static GretlToolItem series_items[] = {
    { N_("Add..."), GTK_STOCK_ADD,    G_CALLBACK(sheet_show_popup),
      SHEET_ADD_BTN },
    { N_("Apply"),  GTK_STOCK_APPLY,  G_CALLBACK(get_data_from_sheet),
      SHEET_APPLY_BTN }
};

static GretlToolItem scalar_items[] = {
    { N_("Add..."),  GTK_STOCK_ADD,    G_CALLBACK(add_scalar_callback), 0 }
};

static int n_series_items = G_N_ELEMENTS(series_items);
static int n_scalar_items = G_N_ELEMENTS(scalar_items);

const gchar *matrix_ui =
    "<ui>"
    "  <menubar>"
    "    <menu action='View'>"
    "      <menuitem action='ViewFmt'/>"
    "      <menuitem action='ViewProps'/>"
    "    </menu>"
    "    <menu action='Fill'>"
    "      <menuitem action='FillIdentity'/>"
    "      <menuitem action='FillUniform'/>"
    "      <menuitem action='FillNormal'/>"
    "    </menu>"
    "    <menu action='Transform'>"
    "      <menuitem action='XTX'/>"
    "      <menuitem action='Transpose'/>"
    "      <menuitem action='Cholesky'/>"
    "      <menuitem action='Invert'/>"
    "      <menuitem action='ScalarMult'/>"
    "      <menuitem action='ScalarDiv'/>"
    "    </menu>"
    "  </menubar>"
    "</ui>";

static GtkActionEntry matrix_items[] = {
    { "View", NULL, N_("_View"), NULL, NULL, NULL },
    { "ViewFmt", NULL, N_("_Format..."), NULL, NULL, G_CALLBACK(matrix_format_callback) },
    { "ViewProps", NULL, N_("_Properties"), NULL, NULL, G_CALLBACK(matrix_props_callback) },
    { "Fill", NULL, N_("_Fill"), NULL, NULL, NULL },
    { "FillIdentity", NULL, N_("_Identity matrix"), NULL, NULL, G_CALLBACK(matrix_fill_callback) },
    { "FillUniform", NULL, N_("_Uniform random"), NULL, NULL, G_CALLBACK(matrix_fill_callback) },
    { "FillNormal", NULL, N_("_Normal random"), NULL, NULL, G_CALLBACK(matrix_fill_callback) },
    { "Transform", NULL, N_("_Transform"), NULL, NULL, NULL },
    { "XTX", NULL, N_("_X'X"), NULL, NULL, G_CALLBACK(matrix_edit_callback) },
    { "Transpose", NULL, N_("_Transpose"), NULL, NULL, G_CALLBACK(matrix_edit_callback) },
    { "Cholesky", NULL, N_("_Cholesky"), NULL, NULL, G_CALLBACK(matrix_edit_callback) },
    { "Invert", NULL, N_("_Invert"), NULL, NULL, G_CALLBACK(matrix_edit_callback) },
    { "ScalarMult", NULL, N_("_Multiply by scalar"), NULL, NULL, G_CALLBACK(matrix_edit_callback) },
    { "ScalarDiv", NULL, N_("_Divide by scalar"), NULL, NULL, G_CALLBACK(matrix_edit_callback) }
};

static Spreadsheet *scalars_sheet;

#define sheet_is_modified(s) (s->flags & SHEET_MODIFIED)

#define editing_scalars(s) (s->cmd == SHEET_EDIT_SCALARS)

static void sheet_set_modified (Spreadsheet *sheet, gboolean s)
{
    if (s) {
	sheet->flags |= SHEET_MODIFIED;
    } else {
	sheet->flags &= ~SHEET_MODIFIED;
	free(sheet->edits);
	sheet->edits = NULL;
	free(sheet->inserts);
	sheet->inserts = NULL;
    }

    if (sheet->save != NULL) {
	gtk_widget_set_sensitive(sheet->save, s);
    } else if (sheet->apply != NULL) {
	gtk_widget_set_sensitive(sheet->apply, s);
    }
}

/* record the fact that the value at @col and @row as been
   changed */

static int sheet_push_edit (Spreadsheet *sheet, int col, int row)
{
    int n_edits = (sheet->edits != NULL)? sheet->edits[0] : 0;
    int *ed = realloc(sheet->edits, (n_edits + 3) * sizeof *ed);

    if (ed == NULL) {
	nomem();
	return E_ALLOC;
    } else {
	sheet->edits = ed;
	n_edits += 2;
	sheet->edits[0] = n_edits;
	sheet->edits[n_edits-1] = col;
	sheet->edits[n_edits] = row;
	sheet_set_modified(sheet, TRUE);
	return 0;
    }
}

/* record the fact that a new row (observation) has been inserted
   at @row */

static int sheet_push_insert (Spreadsheet *sheet, int row)
{
    int *test = gretl_list_append_term(&sheet->inserts, row);

    if (test == NULL) {
	nomem();
	sheet->inserts = NULL;
	return E_ALLOC;
    } else {
	sheet_set_modified(sheet, TRUE);
	return 0;
    }
}

static void sheet_revise_edit_rows (Spreadsheet *sheet, int new_row)
{
    int n_edits = (sheet->edits != NULL)? sheet->edits[0] : 0;

    if (n_edits > 0) {
	int i;

	for (i=2; i<=n_edits; i+=2) {
	    if (sheet->edits[i] >= new_row) {
		sheet->edits[i] += 1;
	    }
	}
    }
}

static int data_column_edited (Spreadsheet *sheet, int col)
{
    if (sheet->edits != NULL) {
	int j;

	for (j=1; j<sheet->edits[0]; j+=2) {
	    if (sheet->edits[j] == col) {
		return 1;
	    }
	}
    }

    return 0;
}

static int data_cell_edited (Spreadsheet *sheet, int col, int row)
{
    if (sheet->edits != NULL) {
	int j;

	for (j=1; j<sheet->edits[0]; j+=2) {
	    if (sheet->edits[j] == col && sheet->edits[j+1] == row) {
		return 1;
	    }
	}
    }

    return 0;
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
    if (sheet->matrix != NULL) {
	gchar *pstr = gtk_tree_path_to_string(path);
	gint cnum = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(column), "colnum"));

	sprintf(sheet->location, "%d, %d", atoi(pstr) + 1, cnum);
	g_free(pstr);
    } else {
	GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(sheet->view));
	const gchar *cstr = gtk_tree_view_column_get_title(column);
	GtkTreeIter iter;
	char tmp[VNAMELEN];
	gchar *rstr;

	single_underscores(tmp, cstr);
	gtk_tree_model_get_iter(model, &iter, path);
	gtk_tree_model_get(model, &iter, 0, &rstr, -1);
	sprintf(sheet->location, "%s, %s", tmp, rstr);
	g_free(rstr);
    }

    gtk_statusbar_pop(GTK_STATUSBAR(sheet->locator), sheet->cid);
    gtk_statusbar_push(GTK_STATUSBAR(sheet->locator),
		       sheet->cid, sheet->location);
}

static void set_treeview_column_number (GtkTreeViewColumn *col, int j)
{
    g_object_set_data(G_OBJECT(col), "colnum", GINT_TO_POINTER(j));
}

static int get_treeview_column_number (GtkTreeViewColumn *col)
{
    if (col != NULL) {
	gpointer p = g_object_get_data(G_OBJECT(col), "colnum");

	if (p != NULL) {
	    return GPOINTER_TO_INT(p);
	}
    }

    return 0;
}

static GtkCellRenderer *get_sheet_renderer (Spreadsheet *sheet,
					    int colnum)
{
    return colnum == 0 ? sheet->textcell : sheet->datacell;
}

static void move_to_next_column (Spreadsheet *sheet, GtkTreePath *path,
				 GtkTreeViewColumn *col)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    int colnum = get_treeview_column_number(col);
    GtkTreeViewColumn *newcol;
    int colmin = 1;

#if CELLDEBUG
    fprintf(stderr, "move_to_next_column: sheet->next = %d\n", sheet->next);
#endif

    if (sheet->next == NEXT_LEFT) {
	colnum--;
    } else {
	colnum++;
    }

#if CELLDEBUG
    fprintf(stderr, "move_to_next_column: colnum = %d\n", colnum);
#endif

    if (editing_scalars(sheet)) {
	colmin = 0;
    }

    if (colnum < colmin || colnum > sheet->datacols) {
	return;
    }

    newcol = gtk_tree_view_get_column(view, colnum);

    if (newcol != NULL) {
	GtkCellRenderer *r = get_sheet_renderer(sheet, colnum);

	gtk_tree_view_set_cursor_on_cell(view, path, newcol, r, FALSE);
    }
}

static void move_to_next_cell (Spreadsheet *sheet, GtkTreePath *path,
			       GtkTreeViewColumn *column)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    gint newrow;

    if (sheet->next == NEXT_SAME) {
	sheet->next = NEXT_DOWN;
	return;
    }

    if (sheet->next == NEXT_RIGHT) {
#if CELLDEBUG
	fprintf(stderr, "move_to_next_cell: got NEXT_RIGHT\n");
#endif
	move_to_next_column(sheet, path, column);
	return;
    }

    if (sheet->next == NEXT_UP) {
	newrow = gtk_tree_path_get_indices(path)[0] - 1;
    } else {
	newrow = gtk_tree_path_get_indices(path)[0] + 1;
    }

#if CELLDEBUG
    fprintf(stderr, "move_to_next_cell: newrow = %d\n", newrow);
#endif

    if (newrow >= 0 && newrow < sheet->datarows) {
	GtkTreePath *newpath;
	gchar pstr[12];

	sprintf(pstr, "%d", newrow);
	newpath = gtk_tree_path_new_from_string(pstr);
	if (newpath != NULL) {
	    int colnum = get_treeview_column_number(column);
	    GtkCellRenderer *r = get_sheet_renderer(sheet, colnum);

	    gtk_tree_view_set_cursor_on_cell(view, newpath, column,
					     r, FALSE);
	    gtk_tree_path_free(newpath);
	}
    } else {
#if CELLDEBUG
	fprintf(stderr, "move_to_next_cell: calling move_to_next_column\n");
#endif
	move_to_next_column(sheet, path, column);
    }

    sheet->next = NEXT_DOWN;
}

static void update_sheet_matrix (Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkTreeModel *model;
    gchar *numstr;
    double x;
    int i, j;

    model = gtk_tree_view_get_model(view);
    gtk_tree_model_get_iter_first(model, &iter);

    for (i=0; i<sheet->datarows; i++) {
	for (j=0; j<sheet->datacols; j++) {
	    gtk_tree_model_get(model, &iter, j + 1, &numstr, -1);
	    if (!strcmp(numstr, "nan")) {
		x = NADBL;
	    } else {
		x = atof(numstr);
	    }
	    gretl_matrix_set(sheet->matrix, i, j, x);
	    g_free(numstr);
	}
	gtk_tree_model_iter_next(model, &iter);
    }
}

static void maybe_update_column_names (Spreadsheet *sheet)
{
    if (sheet->flags & SHEET_COLNAME_MOD) {
	gretl_matrix *M = sheet->oldmat;
	char **cnames;

	cnames = strings_array_new_with_length(M->cols, 13);

	if (cnames == NULL) {
	    nomem();
	} else {
	    GtkTreeViewColumn *col;
	    const char *title;
	    int j, err = 0;

	    for (j=0; j<M->cols && !err; j++) {
		col = gtk_tree_view_get_column(GTK_TREE_VIEW(sheet->view),
					       j + 1);
		if (col != NULL) {
		    title = gtk_tree_view_column_get_title(col);
		    if (title != NULL) {
			single_underscores(cnames[j], title);
		    } else {
			err = 1;
		    }
		}
	    }
	    if (!err) {
		gretl_matrix_set_colnames(M, cnames);
	    }
	}

	sheet->flags &= ~SHEET_COLNAME_MOD;
    }
}

/* in case we're editing a pre-existing matrix, carry the
   modifications back */

static void update_saved_matrix (Spreadsheet *sheet)
{
    if (sheet->oldmat != NULL) {
	gretl_matrix *targ = sheet->oldmat;
	gretl_matrix *src = sheet->matrix;

	if (targ->rows == src->rows && targ->cols == src->cols) {
	    gretl_matrix_copy_values(targ, src);
	} else if (targ->rows * targ->cols == src->rows * src->cols) {
	    gretl_matrix_destroy_info(targ);
	    targ->rows = src->rows;
	    targ->cols = src->cols;
	    gretl_matrix_copy_values(targ, src);
	} else {
	    gretl_matrix *m = gretl_matrix_copy(src);

	    if (m == NULL) {
		nomem();
	    } else {
		user_matrix_replace_matrix_by_name(sheet->mname, m);
		sheet->oldmat = m;
	    }
	}

	maybe_update_column_names(sheet);

	/* record the fact that a matrix has been changed */
	mark_session_changed();
    }
}

static void update_sheet_matrix_element (Spreadsheet *sheet,
					 const gchar *new_text,
					 const gchar *path_string,
					 int colnum)
{
    int i = atoi(path_string);
    int j = colnum - 1;
    double x = atof(new_text);

    if (!strcmp(new_text, "nan")) {
	x = NADBL;
    } else {
	x = atof(new_text);
    }

    gretl_matrix_set(sheet->matrix, i, j, x);
}

static void
maybe_update_store (Spreadsheet *sheet, const gchar *new_text,
		    const gchar *path_string)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeViewColumn *col = NULL;
    GtkTreePath *path;
    GtkTreeIter iter;
    gchar *old_text;
    gint colnum;

    gtk_tree_view_get_cursor(view, NULL, &col);
    if (col == NULL) {
	return;
    }

    colnum = get_treeview_column_number(col);
    path = gtk_tree_path_new_from_string(path_string);
    if (!gtk_tree_model_get_iter(model, &iter, path)) {
	fprintf(stderr, "get_iter failed for path '%s'!\n", path_string);
	return;
    }
    gtk_tree_model_get(model, &iter, colnum, &old_text, -1);

#if CELLDEBUG
    fprintf(stderr, "maybe_update_store: pathstr=%s, colnum=%d, old='%s', new='%s'\n",
	    path_string, colnum, old_text, new_text);
#endif

    if (old_text == NULL || strcmp(old_text, new_text)) {
	gtk_list_store_set(GTK_LIST_STORE(model), &iter,
			   colnum, new_text, -1);

	if (sheet->matrix != NULL) {
	    update_sheet_matrix_element(sheet, new_text, path_string, colnum);
	    sheet_set_modified(sheet, TRUE);
	} else if (editing_scalars(sheet)) {
	    update_scalars_from_sheet(sheet);
	} else {
	    sheet_push_edit(sheet, colnum, atoi(path_string));
	}
    }

    move_to_next_cell(sheet, path, col);
    gtk_tree_path_free(path);
    g_free(old_text);
}

static int scalar_try_genr (gchar **ps)
{
    double x;
    int err = 0;

    x = generate_scalar(*ps + 1, dataset, &err);

    if (!err) {
	g_free(*ps);
	*ps = g_strdup_printf("%.*g", DBL_DIG, x);
    }

    return err;
}

/* right now, this is only used for renaming or adding a scalar */

static void sheet_text_cell_edited (GtkCellRendererText *cell,
				    const gchar *path_string,
				    const gchar *user_text,
				    Spreadsheet *sheet)
{
    int err = 0;

    if (*user_text != '\0') {
	err = gui_validate_varname_strict(user_text,
					  GRETL_TYPE_DOUBLE,
					  sheet->win);
    }

    if (err) {
	GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
	GtkTreePath *path;
	GtkTreeViewColumn *column;

	path = gtk_tree_path_new_from_string(path_string);
	column = gtk_tree_view_get_column(view, 0);
	gtk_tree_view_set_cursor_on_cell(view, path, column,
					 (GtkCellRenderer *) cell,
					 TRUE);
	if (sheet->entry != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(sheet->entry), user_text);
	    gtk_editable_select_region(GTK_EDITABLE(sheet->entry), 0, -1);
	}
	gtk_tree_path_free(path);
    } else {
	maybe_update_store(sheet, user_text, path_string);
    }
}

static void sheet_cell_edited (GtkCellRendererText *cell,
			       const gchar *path_string,
			       const gchar *user_text,
			       Spreadsheet *sheet)
{
    gchar *new_text = NULL;
    int err = 0;

#if CELLDEBUG
    fprintf(stderr, "*** sheet_cell_edited\n");
#endif

    gretl_error_clear();

    if (!strcmp(user_text, "na") || !strcmp(user_text, "NA")) {
	/* allow conversion to missing or NaN */
	if (sheet->matrix != NULL) {
	    new_text = g_strdup("nan");
	} else {
	    new_text = g_strdup("");
	}
    } else {
	new_text = g_strdup(user_text);
	if (*new_text == '=' && editing_scalars(sheet)) {
	    err = scalar_try_genr(&new_text);
	} else {
	    if (sheet->flags & SHEET_USE_COMMA) {
		/* accept point also: convert to locale */
		gretl_charsub(new_text, '.', ',');
	    }
	    err = check_atof(new_text);
	}
	if (err) {
	    errbox(gretl_errmsg_get());
	    g_free(new_text);
	}
    }

    if (!err && new_text != NULL) {
	maybe_update_store(sheet, new_text, path_string);
	g_free(new_text);
    }
}

static void cell_edited_callback (GtkTreeViewColumn *column,
				  const gchar *path_string,
				  const gchar *user_text,
				  Spreadsheet *sheet)

{
    if (editing_scalars(sheet) && get_treeview_column_number(column) == 0) {
	sheet_text_cell_edited(NULL, path_string, user_text, sheet);
    } else {
	sheet_cell_edited(NULL, path_string, user_text, sheet);
    }
}

static void
spreadsheet_scroll_to_new_col (Spreadsheet *sheet, GtkTreeViewColumn *column)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreePath *path;
    gchar *pstr;
    GtkAdjustment *adj;
    GtkWidget *sw;

    pstr = g_strdup("0");
    path = gtk_tree_path_new_from_string(pstr);
    gtk_tree_view_set_cursor(view, path, column, TRUE); /* ?? */
    sw = gtk_widget_get_ancestor(GTK_WIDGET(view), GTK_TYPE_BIN);
    if (sw != NULL) {
	adj = gtk_scrolled_window_get_hadjustment(GTK_SCROLLED_WINDOW(sw));
	gtk_adjustment_set_value(adj, gtk_adjustment_get_upper(adj));
    }
    gtk_tree_path_free(path);
    g_free(pstr);
}

static GtkTreeViewColumn *
add_treeview_column_with_title (Spreadsheet *sheet, const char *name)
{
    GtkTreeViewColumn *column;
    gint cols, colnum;

    column = gtk_tree_view_column_new();
    gtk_tree_view_column_set_title(column, name);
    set_up_sheet_column(column, get_data_col_width(), TRUE);

    cols = gtk_tree_view_insert_column(GTK_TREE_VIEW(sheet->view), column, -1);
    colnum = cols - 1;

    gtk_tree_view_column_pack_start(column, sheet->datacell, TRUE);
    gtk_tree_view_column_set_attributes(column,
					sheet->datacell,
					"text", colnum,
					NULL);
    set_treeview_column_number(column, colnum);

    return column;
}

static int real_add_new_series (Spreadsheet *sheet, const char *varname)
{
    GtkTreeViewColumn *column;
    gchar *tmp = NULL;

    if (add_data_column(sheet)) {
	return 1;
    }

#if SSDEBUG
    fprintf(stderr, "real_add_new_series, after add_data_column: sheet->totcols=%d\n",
	    sheet->totcols);
#endif

    tmp = double_underscores_new(varname);
    column = add_treeview_column_with_title(sheet, tmp);
    g_free(tmp);

    /* scroll to editing position if need be */
    spreadsheet_scroll_to_new_col(sheet, column);

    sheet_set_modified(sheet, TRUE);

    return 0;
}

static void add_scalar_callback (GtkWidget *w, Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeViewColumn *column;
    GtkTreePath *path;
    gchar *pstr;
    GtkListStore *store;
    GtkTreeIter iter;
    GtkAdjustment *adj;
    GtkWidget *sw;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_list_store_append(store, &iter);

    gtk_list_store_set(store, &iter, 0, "", 1, "", 2, sheet->pbuf, -1);
    sheet->datarows += 1;

    pstr = g_strdup_printf("%d", sheet->datarows - 1);
    path = gtk_tree_path_new_from_string(pstr);
    column = gtk_tree_view_get_column(view, 0);
    gtk_tree_view_set_cursor(view, path, column, TRUE);
    gtk_tree_path_free(path);
    g_free(pstr);

    sw = gtk_widget_get_ancestor(GTK_WIDGET(view), GTK_TYPE_BIN);
    if (sw != NULL) {
	adj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(sw));
	gtk_adjustment_set_value(adj, gtk_adjustment_get_upper(adj));
    }
}

static void
spreadsheet_scroll_to_foot (Spreadsheet *sheet, int row, int col)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreePath *path;
    gchar *pstr;
    GtkTreeViewColumn *column;
    GtkAdjustment *adj;
    GtkWidget *sw;

    pstr = g_strdup_printf("%d", row);
    path = gtk_tree_path_new_from_string(pstr);
    column = gtk_tree_view_get_column(view, col);
    gtk_tree_view_set_cursor(view, path, column, FALSE);
    sw = gtk_widget_get_ancestor(GTK_WIDGET(view), GTK_TYPE_BIN);
    if (sw != NULL) {
	adj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(sw));
	gtk_adjustment_set_value(adj, gtk_adjustment_get_upper(adj));
    }
    gtk_tree_path_free(path);
    g_free(pstr);
}

static void
real_add_new_obs (Spreadsheet *sheet, const char *obsname, int n)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    gint rownum = 0;
    gint oldrows = sheet->datarows;
    GtkListStore *store;
    gchar *pstr = NULL;
    GtkTreeIter iter;
    gchar rowlabel[10];
    gint i, j;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));

    if (sheet->point == SHEET_AT_END) {
	rownum = sheet->datarows - 1;
	pstr = g_strdup_printf("%d", rownum);
	for (i=0; i<n; i++) {
	    gtk_list_store_append(store, &iter);
	}
    } else if (sheet->point == SHEET_AT_POINT) {
	GtkTreePath *path;

	gtk_tree_view_get_cursor(view, &path, NULL);
	if (path != NULL) {
	    rownum = gtk_tree_path_get_indices(path)[0];
	    gtk_tree_model_get_iter(GTK_TREE_MODEL(store), &iter, path);
	    gtk_list_store_insert(store, &iter, rownum);
	    gtk_tree_path_free(path);
	}
    } else {
	return;
    }

    if (dataset->markers && obsname != NULL) {
	gtk_list_store_set(store, &iter, 0, obsname, -1);
    } else if (sheet->point == SHEET_AT_END) {
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter,
					    pstr);
	for (j=0; j<n; j++) {
	    gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	    get_obs_string(rowlabel, sheet->datarows + j, dataset);
	    gtk_list_store_set(store, &iter, 0, rowlabel, -1);
	}
    }

    sheet->datarows += n;

    if (pstr != NULL) {
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter,
					    pstr);
	for (j=0; j<n; j++) {
	    gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	    for (i=1; i<=sheet->datacols; i++) {
		gtk_list_store_set(store, &iter, i, "", -1);
	    }
	}
    } else {
	for (i=1; i<=sheet->datacols; i++) {
	    gtk_list_store_set(store, &iter, i, "", -1);
	    sheet_push_edit(sheet, i, rownum);
	}
    }

    if (sheet->point == SHEET_AT_POINT) {
	sheet_revise_edit_rows(sheet, rownum);
	if (!dataset->markers) {
	    for (i=rownum; i<sheet->datarows; i++) {
		get_obs_string(rowlabel, i, dataset);
		gtk_list_store_set(store, &iter, 0, rowlabel, -1);
		gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	    }
	}
    }

    if (sheet->point == SHEET_AT_END) {
	spreadsheet_scroll_to_foot(sheet, oldrows, 1);
    } else {
	GtkTreePath *path;
	GtkTreeIter insiter;

	pstr = g_strdup_printf("%d", rownum);
	path = gtk_tree_path_new_from_string(pstr);
	gtk_tree_model_get_iter(GTK_TREE_MODEL(store), &insiter, path);
	gtk_list_store_set(store, &insiter, 1, "", -1);
	gtk_tree_path_free(path);
	sheet_push_insert(sheet, rownum);
    }

    if (pstr != NULL) {
	g_free(pstr);
    }

    sheet_set_modified(sheet, TRUE);
}

static void name_new_var (GtkWidget *widget, dialog_t *dlg)
{
    Spreadsheet *sheet = (Spreadsheet *) edit_dialog_get_data(dlg);
    GtkWidget *parent = edit_dialog_get_window(dlg);
    const gchar *buf;
    char varname[VNAMELEN];

    buf = edit_dialog_get_text(dlg);

    if (buf == NULL || gui_validate_varname(buf,
					    GRETL_TYPE_SERIES,
					    parent)) {
	return;
    }

    *varname = 0;
    strncat(varname, buf, VNAMELEN - 1);

    edit_dialog_close(dlg);

    if (real_add_new_series(sheet, varname)) {
	nomem();
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

    edit_dialog_close(dlg);
    real_add_new_obs(sheet, obsmarker, 1);
}

static void name_var_dialog (Spreadsheet *sheet)
{
    gchar *msg;

    msg = g_strdup_printf(_("Enter name for new variable\n"
			    "(max. %d characters)"),
			  VNAMELEN - 1);

    edit_dialog(0, _("gretl: name variable"),
		msg, NULL, name_new_var, sheet,
		VARCLICK_NONE, sheet->win);
    g_free(msg);
}

static void new_case_dialog (Spreadsheet *sheet)
{
    gchar *msg;

    msg = g_strdup_printf(_("Enter case marker for new obs\n"
			    "(max. %d characters)"),
			  OBSLEN - 1);

    edit_dialog(0, _("gretl: case marker"),
		msg, NULL, name_new_obs, sheet,
		VARCLICK_NONE, sheet->win);
    g_free(msg);
}

static void name_matrix_col (GtkWidget *widget, dialog_t *dlg)
{
    GtkTreeViewColumn *col = (GtkTreeViewColumn *) edit_dialog_get_data(dlg);
    GtkWidget *parent = edit_dialog_get_window(dlg);
    const gchar *buf, *old;
    char tmp[13], colname[24];

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL || gui_validate_varname(buf,
					    GRETL_TYPE_NONE,
					    parent)) {
	return;
    }

    *tmp = '\0';
    strncat(tmp, buf, 12);
    double_underscores(colname, tmp);

    edit_dialog_close(dlg);

    old = gtk_tree_view_column_get_title(col);

    if (strcmp(old, colname)) {
	Spreadsheet *sheet = g_object_get_data(G_OBJECT(col), "sheet");

	gtk_tree_view_column_set_title(col, colname);
	sheet->flags |= SHEET_COLNAME_MOD;
	sheet_set_modified(sheet, TRUE);
    }
}

static void name_column_dialog (GtkTreeViewColumn *column,
				Spreadsheet *sheet)
{
    edit_dialog(0, _("gretl: name column"),
		_("Enter name for column\n"
		  "(max. 12 characters)"),
		gtk_tree_view_column_get_title(column),
		name_matrix_col, column,
		VARCLICK_NONE, sheet->win);
}

static GtkListStore *make_sheet_liststore (Spreadsheet *sheet)
{
    GtkListStore *store;
    GType *types;
    int ncols;
    int i;

    /* obs col, data cols */
    ncols = sheet->totcols = sheet->datacols + 1;

    if (editing_scalars(sheet)) {
	/* allow for trash-can column */
	ncols++;
    }

    types = mymalloc(ncols * sizeof *types);
    if (types == NULL) {
	return NULL;
    }

    for (i=0; i<sheet->totcols; i++) {
	types[i] = G_TYPE_STRING;
    }

    if (editing_scalars(sheet)) {
	types[i] = GDK_TYPE_PIXBUF;
    }

    store = gtk_list_store_newv(ncols, types);
    free(types);

    return store;
}

/* This is relatively complex because, so far as I can tell, you can't
   append or insert additional columns in a GtkListStore: we have to
   create a whole new liststore and copy the old info across.
*/

static int add_data_column (Spreadsheet *sheet)
{
    GtkListStore *old_store, *new_store;
    GtkTreeIter old_iter, new_iter;
    gint row;

    sheet->datacols += 1;

    /* get pointers to original and new stores */
    new_store = make_sheet_liststore(sheet);
    if (new_store == NULL) {
	sheet->datacols -= 1;
	sheet->totcols -= 1;
	return E_ALLOC;
    }

    old_store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sheet->view)));

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

static void sheet_get_scalar (GtkWidget *w, dialog_t *dlg)
{
    double x, *px = (double *) edit_dialog_get_data(dlg);
    const gchar *buf = edit_dialog_get_text(dlg);
    int err = 0;

    if (buf != NULL) {
	x = gui_double_from_string(buf, &err);
    }

    if (buf == NULL || err) {
	edit_dialog_reset(dlg);
    } else {
	*px = x;
	edit_dialog_close(dlg);
    }
}

static void matrix_edit_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    Spreadsheet *sheet = (Spreadsheet *) data;
    double x = NADBL;
    int err = 0;

    if (!strcmp(s, "ScalarMult") || !strcmp(s, "ScalarDiv")) {
	int cancel = 0;

	blocking_edit_dialog(0, _("gretl: specify scalar"),
			     _("Enter a numerical value"),
			     NULL, sheet_get_scalar, &x,
			     VARCLICK_NONE, sheet->win,
			     &cancel);
	if (cancel || na(x)) {
	    return;
	}
    }

    if (!strcmp(s, "XTX")) {
	err = matrix_XTX_in_place(sheet->matrix);
    } else if (!strcmp(s, "Transpose")) {
	err = matrix_transpose_in_place(sheet->matrix);
    } else if (!strcmp(s, "Cholesky")) {
	err = matrix_cholesky_in_place(sheet->matrix);
    } else if (!strcmp(s, "Invert")) {
	err = matrix_invert_in_place(sheet->matrix);
    } else if (!strcmp(s, "ScalarMult")) {
	gretl_matrix_multiply_by_scalar(sheet->matrix, x);
    } else if (!strcmp(s, "ScalarDiv")) {
	gretl_matrix_divide_by_scalar(sheet->matrix, x);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	update_sheet_from_matrix(sheet);
    }
}

static void matrix_props_callback (GtkAction *action, gpointer data)
{
    Spreadsheet *sheet = (Spreadsheet *) data;

    view_matrix_properties(sheet->matrix, sheet->mname);
}

static void matrix_fill_callback (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    Spreadsheet *sheet = (Spreadsheet *) data;
    gretl_matrix *A = sheet->matrix;
    int mindim;

    if (!strcmp(s, "FillIdentity")) {
	mindim = (A->rows < A->cols)? A->rows : A->cols;
	gretl_matrix_inscribe_I(A, 0, 0, mindim);
    } else if (!strcmp(s, "FillUniform")) {
	gretl_matrix_random_fill(A, D_UNIFORM);
    } else if (!strcmp(s, "FillNormal")) {
	gretl_matrix_random_fill(A, D_NORMAL);
    }

    update_sheet_from_matrix(sheet);
}

static void sheet_add_obs_direct (Spreadsheet *sheet)
{
    if (dataset->markers) {
	new_case_dialog(sheet);
    } else if (sheet->point == SHEET_AT_END) {
	int n = add_obs_dialog(NULL, 1, OPT_NONE, sheet->win);

	if (n > 0) {
	    real_add_new_obs(sheet, NULL, n);
	}
    } else {
	real_add_new_obs(sheet, NULL, 1);
    }
}

static void popup_sheet_add_obs (GtkWidget *w, Spreadsheet *sheet)
{
    sheet->point = SHEET_AT_END;
    sheet_add_obs_direct(sheet);
}

static void popup_sheet_insert_obs (GtkWidget *w, Spreadsheet *sheet)
{
    sheet->point = SHEET_AT_POINT;
    sheet_add_obs_direct(sheet);
}

static void popup_sheet_add_var (GtkWidget *w, Spreadsheet *sheet)
{
    name_var_dialog(sheet);
}

static void sheet_delete_scalar (Spreadsheet *sheet, GtkTreePath *path)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkListStore *store;
    GtkTreeIter iter;
    gchar *vname;
    int err;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_tree_model_get_iter(GTK_TREE_MODEL(store), &iter, path);
    gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 0, &vname, -1);

    set_scalar_edit_callback(NULL);

    err = user_var_delete_by_name(vname, NULL);
    if (!err) {
	mark_session_changed();
    }

    gtk_list_store_remove(store, &iter);
    sheet->datarows--;

    set_scalar_edit_callback(scalars_changed_callback);
}

static void build_sheet_popup (Spreadsheet *sheet)
{
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

/* this is connected to the "cursor-changed" signal
   on the treeview */

static void update_cell_position (GtkTreeView *view,
				  Spreadsheet *sheet)
{
    GtkTreePath *path = NULL;
    GtkTreeViewColumn *col = NULL;
    static int i0, j0;

#if CELLDEBUG > 1
    fprintf(stderr, "*** cursor-changed\n");
#endif

    gtk_tree_view_get_cursor(view, &path, &col);

    if (path != NULL && col != NULL) {
	int i = gtk_tree_path_get_indices(path)[0];
	int j = get_treeview_column_number(col);

	if (j > 0 && (i != i0 || j != j0)) {
#if CELLDEBUG > 1
	    fprintf(stderr, " now in cell(%d, %d)\n", i, j);
#endif
	    set_locator_label(sheet, path, col);
	    i0 = i;
	    j0 = j;
	} else {
#if CELLDEBUG > 1
	   fprintf(stderr, " still in cell(%d, %d)\n", i0, j0);
#endif
	}
    }

    if (path != NULL) {
	gtk_tree_path_free(path);
    }
}

/* put modified values from the spreadsheet into the attached
   matrix */

static void update_matrix_from_sheet_full (Spreadsheet *sheet)
{
    update_sheet_matrix(sheet);

    /* in case we're editing a pre-existing matrix, carry the
       modifications back */
    if (sheet->oldmat != NULL) {
	update_saved_matrix(sheet);
    }

    sheet_set_modified(sheet, FALSE);
    set_ok_transforms(sheet);
}

/* callback from uservar.c, for use when a scalar is added,
   deleted or changed by means other than the spreadsheet, and the
   scalars spreadsheet is currently displayed */

static void scalars_changed_callback (void)
{
    Spreadsheet *sheet = scalars_sheet;

    if (sheet != NULL) {
	GList *slist, *tail;
	GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
	GtkTreeIter iter;
	GtkListStore *store;
	gchar vname[VNAMELEN];
	gchar val[32];
	user_var *u;
	double x;

	tail = slist = user_var_list_for_type(GRETL_TYPE_DOUBLE);

	store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
	gtk_list_store_clear(store);
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	sheet->datarows = 0;

	while (tail) {
	    u = tail->data;
	    if (user_var_get_level(u) != 0) {
		continue;
	    }
	    strcpy(vname, user_var_get_name(u)); /* underscores? */
	    x = user_var_get_scalar_value(u);
	    if (na(x)) {
		*val = '\0';
	    } else {
		sprintf(val, sheet->numfmt, sheet->digits, x);
	    }
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 0, vname, 1, val,
			       2, sheet->pbuf, -1);
	    sheet->datarows += 1;
	    tail = tail->next;
	}

	g_list_free(slist);
    }
}

static double scalar_from_sheet_string (const char *s)
{
    double x = NADBL;

    if (*s != '\0' && strcmp(s, "NA")) {
	if (strstr(s, "nan") || strstr(s, "NAN") || strstr(s, "NaN")) {
	    x = NADBL;
	} else if (strstr(s, "inf") || strstr(s, "INF")) {
	    x = (*s == '-') ? -1.0/0.0 : 1.0/0.0;
	} else {
	    x = atof(s);
	}
    }

    return x;
}

static int scalars_differ (double x, double y)
{
    if (isnan(x) && isnan(y)) {
	return 0;
    } else {
	return x != y;
    }
}

/* put modified values from the spreadsheet into the array
   of saved scalars */

static void update_scalars_from_sheet (Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *vname, *val;
    double x;
    int i, err = 0;

    sheet_set_modified(sheet, FALSE);
    set_scalar_edit_callback(NULL);

    model = gtk_tree_view_get_model(view);
    gtk_tree_model_get_iter_first(model, &iter);

    /* deal with added or modified scalars first */
    for (i=0; i<sheet->datarows && !err; i++) {
	gtk_tree_model_get(model, &iter, 0, &vname, 1, &val, -1);
	if (vname != NULL && *vname != '\0') {
	    x = scalar_from_sheet_string(val);
	    if (gretl_is_scalar(vname)) {
		if (scalars_differ(x, gretl_scalar_get_value(vname, NULL))) {
		    gretl_scalar_set_value(vname, x);
		    sheet_set_modified(sheet, TRUE);
		}
	    } else {
		err = gretl_scalar_add(vname, x);
		if (!err) {
		    sheet_set_modified(sheet, TRUE);
		}
	    }
	}
	g_free(vname);
	g_free(val);
	gtk_tree_model_iter_next(model, &iter);
    }

    if (sheet_is_modified(sheet)) {
	mark_session_changed();
	sheet_set_modified(sheet, FALSE);
    }

    set_scalar_edit_callback(scalars_changed_callback);

    if (err) {
	gui_errmsg(err);
    }
}

/* handle the case where observations have been inserted (not just
   appended); we need to shift the original dataset observations
   "down" by one place for each insertion
*/

static void process_inserts_for_series (Spreadsheet *sheet, int v)
{
    double *src;
    int n = sheet->orig_nobs;
    int i, t;

    for (i=1; i<=sheet->inserts[0]; i++) {
	/* the starting point for the values to be moved */
	t = dataset->t1 + sheet->inserts[i];
	src = dataset->Z[v] + t;
	memmove(src + 1, src, (n - t) * sizeof *src);
	n++;
    }
}

/* pull modified values from the data-editing spreadsheet
   into the main dataset */

static void update_dataset_from_sheet (Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkTreeViewColumn *column;
    GtkTreeModel *model;
    int oldv = dataset->v;
    int newvars = sheet->added_vars;
    int newobs = sheet->datarows - sheet->orig_nobs;
    int missobs = 0;
    int i, j, t;

#if SSDEBUG
    fprintf(stderr, "update_dataset_from_sheet: oldv=%d, newvars=%d, newobs=%d\n",
	    oldv, newvars, newobs);
#endif

    model = gtk_tree_view_get_model(view);

    /* first extend the series length, if needed */

    if (newobs > 0) {
	if (dataset_add_observations(dataset, newobs, OPT_A) ||
	    dataset_destroy_hidden_variables(dataset, 0)) {
	    nomem();
	    return;
	}
    }

    /* then add any new series to the dataset */

    if (newvars > 0) {
	int vj, v0 = sheet->varlist[0];
	const gchar *title;

	if (dataset_add_series(dataset, newvars)) {
	    nomem();
	    return;
	}

	for (j=0; j<newvars; j++) {
	    char newname[VNAMELEN];

	    title = NULL;
	    vj = oldv + j;
	    gretl_list_append_term(&sheet->varlist, vj);
	    if (sheet->varlist == NULL) {
		nomem();
		return;
	    }
	    column = gtk_tree_view_get_column(view, v0 + 1 + j);
	    title = gtk_tree_view_column_get_title(column);
	    single_underscores(newname, title);
	    strcpy(dataset->varname[vj], newname);
	    for (i=0; i<dataset->n; i++) {
		dataset->Z[vj][i] = NADBL;
	    }
#if SSDEBUG
	    fprintf(stderr, " added var %d (%s) from column %d\n",
		    vj, newname, v0 + 1 + j);
#endif
	}
    }

    /* now copy data values from spreadsheet, as needed */

    for (j=1; j<=sheet->varlist[0]; j++) {
	int vj = sheet->varlist[j];
	int var_added = (vj >= oldv);

#if SSDEBUG
	fprintf(stderr, " update data for var %d (%s) from column %d?\n",
		vj, dataset->varname[vj], j);
#endif

	if (!var_added && sheet->inserts != NULL) {
	    process_inserts_for_series(sheet, vj);
	}

	if (var_added || data_column_edited(sheet, j)) {
	    gchar *numstr;

#if SSDEBUG
	    fprintf(stderr, " yes, %s\n", var_added ? "added" : "modified");
#endif
	    gtk_tree_model_get_iter_first(model, &iter);
	    for (i=0; i<sheet->datarows; i++) {
		if (var_added || data_cell_edited(sheet, j, i)) {
		    t = dataset->t1 + i;
		    gtk_tree_model_get(model, &iter, j, &numstr, -1);
		    if (*numstr != '\0') {
			dataset->Z[vj][t] = atof(numstr);
		    } else {
			dataset->Z[vj][t] = NADBL;
			missobs = 1;
		    }
		    g_free(numstr);
		}
		gtk_tree_model_iter_next(model, &iter);
	    }
	}
    }

    /* copy observation markers, if relevant */

    if (dataset_has_markers(dataset)) {
	gchar *marker;

	gtk_tree_model_get_iter_first(model, &iter);
	t = dataset->t1;
	for (i=0; i<sheet->datarows; i++) {
	    gtk_tree_model_get(model, &iter, 0, &marker, -1);
	    strcpy(dataset->S[t++], marker);
	    g_free(marker);
	    gtk_tree_model_iter_next(model, &iter);
	}
    }

    register_data(DATA_APPENDED);

    if (missobs) {
	infobox(_("Warning: there were missing observations"));
    }

    sheet_set_modified(sheet, FALSE);
    sheet->added_vars -= newvars; /* record that these are handled */
    sheet->orig_nobs += newobs;   /* and these too */
}

static void matrix_new_name (GtkWidget *w, dialog_t *dlg)
{
    char *newname = (char *) edit_dialog_get_data(dlg);
    GtkWidget *parent = edit_dialog_get_window(dlg);
    const gchar *buf = edit_dialog_get_text(dlg);

    if (buf == NULL || gui_validate_varname(buf,
					    GRETL_TYPE_MATRIX,
					    parent)) {
	edit_dialog_reset(dlg);
    } else {
	*newname = '\0';
	strncat(newname, buf, VNAMELEN - 1);
	edit_dialog_close(dlg);
    }
}

static void matrix_save_as (GtkWidget *w, Spreadsheet *sheet)
{
    char newname[VNAMELEN];
    int cancel = 0;

    blocking_edit_dialog(0, _("gretl: save matrix"),
			 _("Enter a name"), NULL,
			 matrix_new_name, newname,
			 VARCLICK_NONE, sheet->win,
			 &cancel);

    if (!cancel) {
	gretl_matrix *m;
	user_var *u;
	gchar *tmp;

	m = gretl_matrix_copy(sheet->matrix);
	user_var_add_or_replace(newname,
				GRETL_TYPE_MATRIX,
				m);
	strcpy(sheet->mname, newname);

	tmp = g_strdup_printf("gretl: %s", sheet->mname);
	gtk_window_set_title(GTK_WINDOW(sheet->win), tmp);
	g_free(tmp);

	u = get_user_var_by_name(newname);
	g_object_set_data(G_OBJECT(sheet->win), "object", u);

	sheet->oldmat = m;
	maybe_update_column_names(sheet);
	sheet_set_modified(sheet, FALSE);
    }
}

#define new_matrix(s) (s->cmd == SHEET_EDIT_MATRIX && \
                       s->oldmat == NULL)

static void get_data_from_sheet (GtkWidget *w, Spreadsheet *sheet)
{
    if (!sheet_is_modified(sheet) && !new_matrix(sheet)) {
	infobox(_("No changes were made"));
    } else if (sheet->cmd == SHEET_EDIT_MATRIX) {
	update_matrix_from_sheet_full(sheet);
    } else {
	update_dataset_from_sheet(sheet);
    }
}

static void select_first_editable_cell (Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreePath *path;
    GtkTreeViewColumn *column;

    path = gtk_tree_path_new_from_string("0");

    if (editing_scalars(sheet) && n_user_scalars() == 0) {
	column = gtk_tree_view_get_column(view, 0);
	gtk_tree_view_set_cursor_on_cell(view, path, column,
					 sheet->textcell, FALSE);
    } else {
	column = gtk_tree_view_get_column(view, 1);
	gtk_tree_view_set_cursor_on_cell(view, path, column,
					 sheet->datacell, FALSE);
    }

    if (sheet->locator != NULL) {
	set_locator_label(sheet, path, column);
    }

    gtk_tree_path_free(path);
}

static void set_ok_transforms (Spreadsheet *sheet)
{
    const gretl_matrix *m = sheet->matrix;
    int z = gretl_is_zero_matrix(m);
    int s = gretl_matrix_get_structure(m);

    flip(sheet->ui, "/menubar/Fill/FillIdentity", m->rows == m->cols);
    flip(sheet->ui, "/menubar/Transform/ScalarMult", !z);
    flip(sheet->ui, "/menubar/Transform/ScalarDiv", !z);
    flip(sheet->ui, "/menubar/Transform/XTX",
	 s == 0 || (!z && s != GRETL_MATRIX_IDENTITY));
    flip(sheet->ui, "/menubar/Transform/Cholesky",
	 s == GRETL_MATRIX_SYMMETRIC && !z);
    flip(sheet->ui, "/menubar/Transform/Invert", s > 0 && !z &&
	 s != GRETL_MATRIX_IDENTITY);
    flip(sheet->ui, "/menubar/Transform/Transpose",
	 s < GRETL_MATRIX_SYMMETRIC);
}

static void maybe_rename_sheet_cols (Spreadsheet *sheet)
{
    GtkTreeViewColumn *col;
    const char *cname;
    char cstr[16];
    int i;

    for (i=1; i<=sheet->matrix->cols; i++) {
	col = gtk_tree_view_get_column(GTK_TREE_VIEW(sheet->view), i);
	cname = gtk_tree_view_column_get_title(col);
	sprintf(cstr, "%d", i);
	if (strcmp(cname, cstr)) {
	    gtk_tree_view_column_set_title(col, cstr);
	}
    }
}

static int rejig_sheet_cols (Spreadsheet *sheet)
{
    int n = sheet->matrix->cols - sheet->datacols;
    GtkTreeViewColumn *col;
    GtkListStore *store;
    char cstr[16];
    int i;

    sheet->datacols = sheet->matrix->cols;

    store = make_sheet_liststore(sheet);
    if (store == NULL) {
	return E_ALLOC;
    }

    gtk_tree_view_set_model(GTK_TREE_VIEW(sheet->view), GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    if (n > 0) {
	/* add extra columns */
	for (i=0; i<n; i++) {
	    sprintf(cstr, "%d", sheet->datacols - n + i + 1);
	    add_treeview_column_with_title(sheet, cstr);
	}
    } else {
	/* remove surplus columns */
	n = -n;
	for (i=0; i<n; i++) {
	    col = gtk_tree_view_get_column(GTK_TREE_VIEW(sheet->view),
					   sheet->datacols + n - i);
	    gtk_tree_view_remove_column(GTK_TREE_VIEW(sheet->view),
					col);
	}
    }

    for (i=1; i<=sheet->matrix->cols; i++) {
	sprintf(cstr, "%d", i);
	col = gtk_tree_view_get_column(GTK_TREE_VIEW(sheet->view), i);
	gtk_tree_view_column_set_title(col, cstr);
	gtk_tree_view_column_set_clickable(col, TRUE);
	g_object_set_data(G_OBJECT(col), "sheet", sheet);
	g_signal_connect(G_OBJECT(col), "clicked",
			 G_CALLBACK(name_column_dialog), sheet);
    }

    return 0;
}

static int update_sheet_from_matrix (Spreadsheet *sheet)
{
    gchar tmpstr[48];
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkListStore *store;
    GtkTreeIter iter;
    int i, j, resized = 0;
    double x;
    int err = 0;

    if (sheet->matrix->cols != sheet->datacols) {
	err = rejig_sheet_cols(sheet);
	if (err) {
	    return err;
	}
	resized = 1;
    } else {
	maybe_rename_sheet_cols(sheet);
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<sheet->matrix->rows; i++) {
	int i1 = i + 1;

	if (resized) {
	    gtk_list_store_append(store, &iter);
	    sprintf(tmpstr, "%d", i1);
	    gtk_list_store_set(store, &iter, 0, tmpstr, -1);
	} else if (sheet->rownames != NULL) {
	    sprintf(tmpstr, "%d", i1);
	    gtk_list_store_set(store, &iter, 0, tmpstr, -1);
	}

	for (j=0; j<sheet->matrix->cols; j++) {
	    x = gretl_matrix_get(sheet->matrix, i, j);
	    if (x == -0.0) {
		x = 0.0;
	    }
	    sprintf(tmpstr, sheet->numfmt, sheet->digits, x);
	    gtk_list_store_set(store, &iter, j + 1, tmpstr, -1);
	}
	if (i1 < sheet->matrix->rows) {
	    /* there are more rows to handle */
	    if (resized) {
		/* will be handled above on next iteration */
		continue;
	    } else if (i1 >= sheet->datarows) {
		/* need to add another row */
		gtk_list_store_append(store, &iter);
		sprintf(tmpstr, "%d", i1 + 1);
		gtk_list_store_set(store, &iter, 0, tmpstr, -1);
	    } else {
		/* row already present */
		gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	    }
	}
    }

    if (!resized && (sheet->datarows > sheet->matrix->rows)) {
	int rdel = sheet->datarows - sheet->matrix->rows;

	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	for (i=0; i<rdel; i++) {
	    gtk_list_store_remove(store, &iter);
	}
    }

    if (sheet->datarows != sheet->matrix->rows) {
	sheet->datarows = sheet->matrix->rows;
	resized = 1;
    }

    sheet_set_modified(sheet, TRUE);
    set_ok_transforms(sheet);

    if (resized) {
	size_matrix_window(sheet);
    }

    return 0;
}

static int add_matrix_data_to_sheet (Spreadsheet *sheet)
{
    gchar tmpstr[32];
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkListStore *store;
    GtkTreeIter iter;
    double x;
    int i, j;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));

    /* row labels */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    for (i=0; i<sheet->datarows; i++) {
	if (sheet->rownames != NULL) {
	    *tmpstr = '\0';
	    strncat(tmpstr, sheet->rownames[i], 31);
	} else {
	    sprintf(tmpstr, "%d", i+1);
	}
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, tmpstr, -1);
    }

    /* now insert data values */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<sheet->datarows; i++) {
	for (j=0; j<sheet->datacols; j++) {
	    x = gretl_matrix_get(sheet->matrix, i, j);
	    if (isnan(x)) {
		strcpy(tmpstr, "nan");
	    } else {
		sprintf(tmpstr, sheet->numfmt, sheet->digits, x);
	    }
	    gtk_list_store_set(store, &iter, j + 1, tmpstr, -1);
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }

    set_ok_transforms(sheet);

    return 0;
}

static void reformat_sheet_matrix (Spreadsheet *sheet)
{
    gchar tmpstr[32];
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkListStore *store;
    GtkTreeIter iter;
    double x;
    int i, j;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<sheet->datarows; i++) {
	for (j=0; j<sheet->datacols; j++) {
	    x = gretl_matrix_get(sheet->matrix, i, j);
	    if (isnan(x)) {
		strcpy(tmpstr, "nan");
	    } else {
		sprintf(tmpstr, sheet->numfmt, sheet->digits, x);
	    }
	    gtk_list_store_set(store, &iter, j + 1, tmpstr, -1);
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }

    gtk_tree_view_columns_autosize(view);
}

static int add_scalars_to_sheet (Spreadsheet *sheet)
{
    GList *slist, *tail;
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    GtkTreeIter iter;
    GtkListStore *store;
    user_var *u;
    gchar vname[VNAMELEN];
    gchar val[32];
    double x;

    if (sheet->pbuf == NULL) {
	sheet->pbuf = gtk_widget_render_icon(sheet->view, GTK_STOCK_DELETE,
					     GTK_ICON_SIZE_MENU, NULL);
    }

    tail = slist = user_var_list_for_type(GRETL_TYPE_DOUBLE);

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* insert names and values of scalar variables */

    while (tail) {
	u = tail->data;
	strcpy(vname, user_var_get_name(u)); /* underscores? */
	x = user_var_get_scalar_value(u);
	if (na(x)) {
	    if (isnan(x)) {
		strcpy(val, "nan");
	    } else if (x < 0) {
		strcpy(val, "-inf");
	    } else {
		strcpy(val, "inf");
	    }
	} else {
	    sprintf(val, sheet->numfmt, sheet->digits, x);
	}
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, vname, 1, val,
			   2, sheet->pbuf, -1);
	tail = tail->next;
    }

    g_list_free(slist);

    return 0;
}

static void maybe_update_scalars_sheet (Spreadsheet *sheet)
{
    GtkTreeView *view = GTK_TREE_VIEW(sheet->view);
    int ns = n_user_scalars();
    int update = 0;

    if (sheet->datarows != ns) {
	sheet->datarows = ns;
	update = 1;
    } else {
	GList *slist, *tail;
	GtkTreeIter iter;
	GtkTreeModel *model;
	gchar *rowname, *numstr;
	const char *sname;
	user_var *u;
	double sx, rx;

	tail = slist = user_var_list_for_type(GRETL_TYPE_DOUBLE);

	model = gtk_tree_view_get_model(view);
	gtk_tree_model_get_iter_first(model, &iter);

	while (tail) {
	    u = tail->data;
	    gtk_tree_model_get(model, &iter, 0, &rowname, 1, &numstr, -1);
	    rx = atof(numstr);
	    sname = user_var_get_name(u);
	    sx = user_var_get_scalar_value(u);
	    if (strcmp(rowname, sname) || rx != sx) {
		update = 1;
	    }
	    g_free(rowname);
	    g_free(numstr);
	    if (update) {
		break;
	    }
	    gtk_tree_model_iter_next(model, &iter);
	    tail = tail->next;
	}

	g_list_free(slist);
    }

    if (update) {
	GtkListStore *store;

	store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
	gtk_list_store_clear(store);
	add_scalars_to_sheet(sheet);
    }
}

static int add_data_to_sheet (Spreadsheet *sheet, SheetCmd c)
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
    for (t=dataset->t1; t<=dataset->t2; t++) {
	get_obs_string(rowlabel, t, dataset);
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, rowlabel, -1);
    }

    sheet->datarows = dataset->t2 - dataset->t1 + 1;

    /* insert data values */

    if (c == SHEET_NEW_DATASET) {
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	for (t=dataset->t1; t<=dataset->t2; t++) {
	    gtk_list_store_set(store, &iter, 1, "", -1);
	}
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    } else {
	GtkTreeViewColumn *col;
	int vi, numlen, maxwidth;
	char numstr[32];

	for (i=1; i<=sheet->varlist[0]; i++) {
	    maxwidth = 0;
	    vi = sheet->varlist[i];
	    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	    for (t=dataset->t1; t<=dataset->t2; t++) {
		if (na(dataset->Z[vi][t])) {
		    *numstr = '\0';
		} else {
		    sprintf(numstr, sheet->numfmt, sheet->digits,
			    dataset->Z[vi][t]);
		    numlen = strlen(numstr);
		    if (numlen > maxwidth) {
			maxwidth = numlen;
		    }
		}
		gtk_list_store_set(store, &iter, i, numstr, -1);
		gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
	    }
	    col = gtk_tree_view_get_column(view, i);
	    memset(numstr, 0, sizeof numstr);
	    for (t=0; t<maxwidth; t++) {
		numstr[t] = '0';
	    }
	    maxwidth = get_string_width(numstr);
	    gtk_tree_view_column_set_max_width(col, 4 + maxwidth);
	}
    }

    sheet->orig_main_v = dataset->v;

#if SSDEBUG
    fprintf(stderr, " datarows=%d, orig_vars=%d, orig_main_v=%d\n",
	    sheet->datarows, sheet->varlist[0], sheet->orig_main_v);
#endif

    return 0;
}

static gint get_row_label_width (Spreadsheet *sheet)
{
    static gint width;

    if (width == 0) {
	if (sheet != NULL && sheet->rownames != NULL) {
	    const char *s = sheet->rownames[0];
	    int i, len, maxlen = 0;

	    for (i=0; i<sheet->datarows; i++) {
		len = strlen(sheet->rownames[i]);
		if (len > maxlen) {
		    maxlen = len;
		    s = sheet->rownames[i];
		}
	    }
	    width = get_string_width(s);
	} else {
	    width = get_string_width("XXXX");
	}
    }
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

static gint get_name_col_width (void)
{
    static gint width;

    if (width == 0) {
	width = get_string_width("XXXXXXXXXXXXXXX");
    }
    return width;
}

static gint get_data_col_width (void)
{
    static gint width;

    if (width == 0) {
	width = get_string_width("-00.0000000000");
    }
    return width;
}

static gint get_delete_col_width (void)
{
    static gint width;

    if (width == 0) {
	width = get_string_width("XXXXXX");
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

static void commit_and_move_right (Spreadsheet *sheet,
				   const char *s)
{
    GtkTreeViewColumn *col;
    GtkTreePath *path;
    gchar *pathstr;

    gtk_tree_view_get_cursor(GTK_TREE_VIEW(sheet->view), &path, &col);

    if (path == NULL) {
	return;
    }

    pathstr = gtk_tree_path_to_string(path);

#if CELLDEBUG
    fprintf(stderr, "commit_and_move_right: calling sheet_cell_edited\n");
#endif

    cell_edited_callback(col, pathstr, s, sheet);
    g_free(pathstr);
    gtk_tree_path_free(path);
}

static gboolean
catch_sheet_edit_key (GtkWidget *view, GdkEventKey *key,
		      Spreadsheet *sheet)
{
    sheet->next = NEXT_DOWN;

    if (key->keyval == GDK_Up) {
#if CELLDEBUG
	fprintf(stderr, "catch_edit_key: GDK_Up\n");
#endif
	sheet->next = NEXT_UP;
    } else if (key->keyval == GDK_Right && sheet->entry != NULL) {
	int n, pos = gtk_editable_get_position(GTK_EDITABLE(sheet->entry));
	const char *s = gtk_entry_get_text(GTK_ENTRY(sheet->entry));

#if CELLDEBUG
	fprintf(stderr, "catch_edit_key: GDK_Right\n");
#endif

	n = (s != NULL)? strlen(s) : 0;
	if (pos == n) {
	    sheet->next = NEXT_RIGHT;
	    commit_and_move_right(sheet, s);
	    return TRUE;
	}
    } else if (key->keyval == GDK_Return) {
#if CELLDEBUG
	fprintf(stderr, "catch_edit_key: GDK_Return\n");
#endif
    } else if (key->keyval == GDK_Down) {
#if CELLDEBUG
	fprintf(stderr, "catch_edit_key: GDK_Down\n");
#endif
    } else if (key->keyval == GDK_Tab) {
#if CELLDEBUG
	fprintf(stderr, "catch_edit_key: GDK_Tab\n");
#endif
    }

    return FALSE;
}

static void nullify_sheet_entry (gpointer p, Spreadsheet *sheet)
{
#if CELLDEBUG
    fprintf(stderr, "editing entry destroyed\n");
#endif
    sheet->entry = NULL;
}

static void cell_edit_start (GtkCellRenderer *r,
			     GtkCellEditable *ed,
			     gchar *path,
			     Spreadsheet *sheet)
{
#if CELLDEBUG
    fprintf(stderr, "*** editing-started\n");
#endif
    if (GTK_IS_ENTRY(ed)) {
	sheet->entry = GTK_WIDGET(ed);
	g_signal_connect(G_OBJECT(ed), "key-press-event",
			 G_CALLBACK(catch_sheet_edit_key), sheet);
	g_signal_connect(G_OBJECT(ed), "destroy",
			 G_CALLBACK(nullify_sheet_entry), sheet);
    }
}

static void create_sheet_cell_renderers (Spreadsheet *sheet)
{
    GtkCellRenderer *r;

    sheet->textcell = r = gtk_cell_renderer_text_new();

    if (editing_scalars(sheet)) {
	/* editable name */
	g_object_set(r, "ypad", 1,
		     "xalign", 1.0,
		     "editable", TRUE,
		     "mode", GTK_CELL_RENDERER_MODE_EDITABLE,
		     NULL);
	g_signal_connect(r, "editing-started",
			 G_CALLBACK(cell_edit_start), sheet);
	g_signal_connect(r, "edited",
			 G_CALLBACK(sheet_text_cell_edited), sheet);
    } else {
	/* immutable label */
	g_object_set(r, "ypad", 1,
		     "xalign", 1.0,
		     "background", "#DEDEDE",
		     "mode", GTK_CELL_RENDERER_MODE_INERT,
		     NULL);
    }

    sheet->datacell = r = gtk_cell_renderer_text_new();

    g_object_set(r, "ypad", 1,
		 "xalign", 1.0,
		 "family", "Monospace",
		 "editable", TRUE,
		 "mode", GTK_CELL_RENDERER_MODE_EDITABLE,
		 NULL);

    g_signal_connect(r, "editing-started",
		     G_CALLBACK(cell_edit_start), sheet);
    g_signal_connect(r, "edited",
		     G_CALLBACK(sheet_cell_edited), sheet);
}

static void manufacture_keystroke (GtkWidget *widget,
				   GdkEvent *orig,
				   guint uval)
{
    GdkKeymap *keymap = NULL;
    GdkKeymapKey *keys;
    gint n_keys;

#if GTK_MAJOR_VERSION >= 3 || defined(G_OS_WIN32)
    /* with GDK 3, we can't pass NULL for keymap below -- and neither
       (it appears) for GDK 2 on MS Windows
    */
    keymap = gdk_keymap_get_for_display(gdk_display_get_default());
#endif

    if (gdk_keymap_get_entries_for_keyval(keymap, uval, &keys, &n_keys)) {
	guint16 hardware_keycode;
	GdkEvent *event;

	hardware_keycode = keys[0].keycode;
	g_free(keys);

	event = gdk_event_new(GDK_KEY_PRESS);
	event->key.window = g_object_ref(gtk_widget_get_window(widget));
	event->key.hardware_keycode = hardware_keycode;
	event->key.keyval = gdk_unicode_to_keyval(uval);
	event->key.length = 1;
	event->key.send_event = FALSE;
	event->key.time = GDK_CURRENT_TIME;

#if GTK_MAJOR_VERSION >= 3
	/* we now get warning spew if no device is attached */
	gdk_event_set_device(event, gdk_event_get_device(orig));
#endif

	gtk_main_do_event(event);
	gdk_event_free(event);
    }
}

static int numeric_key (guint *kval, Spreadsheet *sheet)
{
    if (*kval >= GDK_KP_0 && *kval <= GDK_KP_9) {
	*kval = GDK_0 + (*kval - GDK_KP_0);
    }

    if (*kval >= GDK_0 && *kval <= GDK_9) {
	return 1;
    } else if (*kval == GDK_minus || *kval == GDK_period) {
	return 1;
    } else if ((sheet->flags & SHEET_USE_COMMA) && *kval == GDK_comma) {
	return 1;
    }

    return 0;
}

#define alpha_key(k) ((k >= GDK_A && k <=  GDK_Z) || (k >= GDK_a && k <= GDK_z))

static gint catch_spreadsheet_key (GtkWidget *view, GdkEvent *event,
				   Spreadsheet *sheet)
{
    guint kval = ((GdkEventKey*) event)->keyval;

#if CELLDEBUG
    fprintf(stderr, "catch_spreadsheet_key: %d (tab=%d)\n", kval, GDK_Tab);
#endif

    if (kval == GDK_Tab) {
	/* FIXME */
	return FALSE;
    }

    if (kval == GDK_Left) {
	if (sheet->cmd != SHEET_EDIT_SCALARS) {
	    GtkTreeViewColumn *col = NULL;

	    gtk_tree_view_get_cursor(GTK_TREE_VIEW(view), NULL, &col);

	    if (col != NULL) {
		if (get_treeview_column_number(col) == 1) {
		    /* if in column 1, don't move left */
		    return TRUE;
		}
	    }

	    return FALSE;
	}
    }

    if (kval == GDK_Up || kval == GDK_Down) {
	GtkTreePath *path = NULL;
	GtkTreeViewColumn *col = NULL;
	int i;

	gtk_tree_view_get_cursor(GTK_TREE_VIEW(view), &path, &col);
	i = (gtk_tree_path_get_indices(path))[0];

	if (kval == GDK_Down && i < sheet->datarows - 1) {
	    gtk_tree_path_next(path);
	} else if (kval == GDK_Up && i > 0) {
	    gtk_tree_path_prev(path);
	}
	gtk_tree_view_set_cursor(GTK_TREE_VIEW(view), path, col,
				 FALSE);
	gtk_tree_path_free(path);
	return TRUE;
    } else if (numeric_key(&kval, sheet)) {
	GtkTreePath *path = NULL;
	GtkTreeViewColumn *col = NULL;

#if CELLDEBUG
	fprintf(stderr, "numeric key: start editing, k = %d\n", kval);
#endif

	gtk_tree_view_get_cursor(GTK_TREE_VIEW(view), &path, &col);
	if (path != NULL && col != NULL) {
	    gtk_tree_view_set_cursor(GTK_TREE_VIEW(view), path, col,
				     TRUE);
	    gtk_tree_path_free(path);
	    manufacture_keystroke(view, event, kval);
	}
    } else if (editing_scalars(sheet) && alpha_key(kval)) {
	GtkTreePath *path = NULL;
	GtkTreeViewColumn *col = NULL;
	gboolean ret = FALSE;

	gtk_tree_view_get_cursor(GTK_TREE_VIEW(view), &path, &col);

	if (path != NULL && col != NULL) {
	    if (get_treeview_column_number(col) >= 1) {
		/* not a text cell */
		ret = TRUE;
	    } else {
		gtk_tree_view_set_cursor(GTK_TREE_VIEW(view), path, col,
					 TRUE);
		manufacture_keystroke(view, event, kval);
	    }
	    gtk_tree_path_free(path);
	}

	return ret;
    }

    return FALSE;
}

static gint catch_spreadsheet_click (GtkWidget *view,
				     GdkEventButton *event,
				     Spreadsheet *sheet)
{
    gint ret = FALSE;

#if CELLDEBUG
    fprintf(stderr, "** catch_spreadsheet_click()\n");
#endif

    if (event->type != GDK_BUTTON_PRESS) {
	return FALSE;
    }

    if (sheet->matrix == NULL && !editing_scalars(sheet) &&
	(right_click(event))) {

	if (sheet->popup == NULL) {
	    build_sheet_popup(sheet);
	}

	gtk_menu_popup(GTK_MENU(sheet->popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
	return TRUE;
    }

    if (event->button == 1) {
	GtkTreePath *path = NULL, *oldpath = NULL;
	GtkTreeViewColumn *column = NULL, *oldcol = NULL;

#if CELLDEBUG
	fprintf(stderr, "Got button 1 click\n");
#endif

	/* where's the cursor at present? */
	gtk_tree_view_get_cursor(GTK_TREE_VIEW(sheet->view),
				 &oldpath, &oldcol);

	if (oldpath != NULL) {
	    if (oldcol != NULL && sheet->entry != NULL) {
		const gchar *txt = gtk_entry_get_text(GTK_ENTRY(sheet->entry));
		gchar *pathstr = gtk_tree_path_to_string(oldpath);

#if CELLDEBUG
		fprintf(stderr, "click: calling sheet_cell_edited\n");
#endif
		cell_edited_callback(oldcol, pathstr, txt, sheet);
		g_free(pathstr);
	    }
	    gtk_tree_path_free(oldpath);
	}

	gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(sheet->view),
				      (gint) event->x,
				      (gint) event->y,
				      &path, &column,
				      NULL, NULL);

	if (path != NULL && column != NULL) {
	    gint colnum = get_treeview_column_number(column);

#if CELLDEBUG
	    fprintf(stderr, "*** Clicked column: colnum = %d\n", colnum);
#endif

	    if (colnum == 0 && !editing_scalars(sheet)) {
		/* don't respond to a click in a non-editable column */
		ret = TRUE;
	    } else if (editing_scalars(sheet) && colnum == 2) {
		/* call to delete scalar */
		sheet_delete_scalar(sheet, path);
	    } else {
		/* activate clicked cell for editing */
		gtk_tree_view_set_cursor(GTK_TREE_VIEW(sheet->view),
					 path, column, TRUE);
		ret = TRUE;
	    }
	}
	gtk_tree_path_free(path);
    }

#if CELLDEBUG
    fprintf(stderr, "catch_spreadsheet_click returning %d\n", ret);
#endif

    return ret;
}

static void sheet_show_popup (GtkWidget *w, Spreadsheet *sheet)
{
    if (sheet->popup == NULL) {
	build_sheet_popup(sheet);
    }

    gtk_menu_popup(GTK_MENU(sheet->popup), NULL, NULL, NULL, NULL,
		   1, gtk_get_current_event_time());
}

static int build_sheet_view (Spreadsheet *sheet)
{
    GtkListStore *store;
    GtkWidget *view;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    gchar *col0str = NULL;
    gchar *tmpstr = NULL;
    gint i, width, colnum;

    if (get_local_decpoint() == ',') {
	sheet->flags |= SHEET_USE_COMMA;
    }

    store = make_sheet_liststore(sheet);
    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), FALSE);
    gtk_tree_view_set_grid_lines(GTK_TREE_VIEW(view), GTK_TREE_VIEW_GRID_LINES_BOTH);

    /* build and attach the (two) cell renderers */
    create_sheet_cell_renderers(sheet);

    /* construct the first (text) column */
    if (editing_scalars(sheet)) {
	col0str = _("Name");
	width = get_name_col_width();
    } else if (sheet->matrix != NULL) {
	width = get_row_label_width(sheet);
    } else {
	width = get_obs_col_width();
    }
    column = gtk_tree_view_column_new_with_attributes(col0str,
						      sheet->textcell,
						      "text", 0,
						      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
    set_up_sheet_column(column, width, FALSE);
    set_treeview_column_number(column, 0);

    width = get_data_col_width();

    if (sheet->matrix != NULL) {
	for (i=1; i<=sheet->datacols; i++) {
	    if (sheet->colnames != NULL) {
		tmpstr = double_underscores_new(sheet->colnames[i-1]);
	    } else {
		tmpstr = g_strdup_printf("%d", i);
	    }
	    column = gtk_tree_view_column_new_with_attributes(tmpstr,
							      sheet->datacell,
							      "text",
							      i,
							      NULL);
	    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	    set_up_sheet_column(column, width, TRUE);
	    set_treeview_column_number(column, i);
	    g_object_set_data(G_OBJECT(column), "sheet", sheet);
	    gtk_tree_view_column_set_clickable(column, TRUE);
	    g_signal_connect(G_OBJECT(column), "clicked",
			     G_CALLBACK(name_column_dialog), sheet);
	    g_free(tmpstr);
	}
	sheet->colnames = NULL;
    } else if (editing_scalars(sheet)) {
	GtkCellRenderer *pixcell;

	/* numerical value */
	column = gtk_tree_view_column_new_with_attributes(_("Value"),
							  sheet->datacell,
							  "text",
							  1,
							  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	set_up_sheet_column(column, width, TRUE);
	set_treeview_column_number(column, 1);
	g_object_set_data(G_OBJECT(column), "sheet", sheet);

	/* deletion icon */
	pixcell = gtk_cell_renderer_pixbuf_new();
	column = gtk_tree_view_column_new_with_attributes(_("Delete"),
							  pixcell,
							  "pixbuf",
							  2,
							  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	gtk_tree_view_column_set_resizable(column, FALSE);
	set_treeview_column_number(column, 2);
	g_object_set_data(G_OBJECT(column), "sheet", sheet);
    } else {
	/* editing data series */
	colnum = 0;
	for (i=1; i<=sheet->varlist[0]; i++) {
	    int vi = sheet->varlist[i];

	    if (series_is_hidden(dataset, vi)) {
		continue;
	    }
	    colnum++;
	    tmpstr = double_underscores_new(dataset->varname[vi]);
	    column = gtk_tree_view_column_new_with_attributes(tmpstr,
							      sheet->datacell,
							      "text",
							      colnum,
							      NULL);
	    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	    set_up_sheet_column(column, width, TRUE);
	    set_treeview_column_number(column, colnum);
	    g_free(tmpstr);
	}
    }

    /* set the selection property on the tree view */
    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_NONE);

    if (sheet->cmd != SHEET_EDIT_SCALARS) {
	g_signal_connect(G_OBJECT(view), "cursor-changed",
			 G_CALLBACK(update_cell_position), sheet);
    }

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

    if (sheet->ui != NULL) {
	g_object_unref(sheet->ui);
    }

    free(sheet->varlist);
    free(sheet->edits);
    free(sheet->inserts);

    if (editing_series(sheet)) {
	set_dataset_locked(FALSE);
    }

    if (sheet->matrix != NULL && sheet->oldmat != NULL) {
	/* delete the copied matrix */
	gretl_matrix_free(sheet->matrix);
    }

    if (editing_scalars(sheet)) {
	scalars_sheet = NULL;
	set_scalar_edit_callback(NULL);
    }

    free(sheet);

    *psheet = NULL;
}

static void free_matrix_sheet (GtkWidget *widget, Spreadsheet *sheet)
{
    if (sheet->popup != NULL) {
	gtk_widget_destroy(sheet->popup);
    }

    free(sheet->varlist);

    if (editing_series(sheet)) {
	set_dataset_locked(FALSE);
    }

    if (sheet->matrix != NULL && sheet->oldmat != NULL) {
	/* delete the copied matrix */
	gretl_matrix_free(sheet->matrix);
    }

    free(sheet);
}

static int sheet_list_empty (Spreadsheet *sheet)
{
    int ret = 0;

    if (sheet->varlist == NULL) {
	free(sheet);
	ret = 1;
    } else if (sheet->varlist[0] == 0) {
	free(sheet->varlist);
	free(sheet);
	ret = 1;
    }

    if (ret) {
	errbox(_("No series to edit"));
    }

    return ret;
}

static Spreadsheet *spreadsheet_new (SheetCmd c, int varnum)
{
    Spreadsheet *sheet;

    sheet = mymalloc(sizeof *sheet);
    if (sheet == NULL) return NULL;

    sheet->view = NULL;
    sheet->win = NULL;
    sheet->locator = NULL;
    sheet->entry = NULL;
    sheet->popup = NULL;
    sheet->ui = NULL;
    sheet->save = NULL;
    sheet->apply = NULL;
    sheet->textcell = NULL;
    sheet->datacell = NULL;
    sheet->pbuf = NULL;
    sheet->datacols = sheet->datarows = 0;
    sheet->totcols = 0;
    sheet->added_vars = 0;
    sheet->orig_main_v = 0;
    sheet->orig_nobs = 0;
    sheet->next = NEXT_DOWN;
    sheet->cid = 0;
    sheet->varlist = NULL;
    sheet->matrix = NULL;
    sheet->oldmat = NULL;
    sheet->colnames = NULL;
    sheet->rownames = NULL;
    sheet->cmd = c;
    sheet->flags = 0;

    sheet->edits = NULL;
    sheet->inserts = NULL;

    sheet->mname[0] = '\0';

    strcpy(sheet->numfmt, "%.*g");

    if (c == SHEET_EDIT_MATRIX || c == SHEET_EDIT_SCALARS) {
	sheet->digits = DBL_DIG;
	return sheet;
    }

    sheet->digits = SERIES_DIGITS;

    sheet->orig_nobs = sample_size(dataset);
    if (sheet->orig_nobs < dataset->n) {
	sheet->flags |= SHEET_SUBSAMPLED;
    }

    if (sheet->cmd == SHEET_NEW_DATASET) {
	sheet->varlist = gretl_list_new(1);
    } else {
	if (sheet->cmd == SHEET_EDIT_VARLIST) {
	    if (varnum > 0) {
		sheet->varlist = gretl_list_new(1);
		if (sheet->varlist != NULL) {
		    sheet->varlist[1] = varnum;
		}
	    } else {
		sheet->varlist = main_window_selection_as_list();
	    }
	} else {
	    sheet->varlist = full_var_list(dataset, NULL);
	}
	if (sheet_list_empty(sheet)) {
	    return NULL;
	}
    }

    if (sheet->varlist == NULL) {
	free(sheet);
	sheet = NULL;
    } else if (sheet->cmd == SHEET_NEW_DATASET) {
	sheet->varlist[1] = 1;
    }

    return sheet;
}

static void empty_dataset_guard (void)
{
    int t, empty = 0, miss = 0;

    if (dataset->v == 2) {
	empty = 1;
	for (t=0; t<dataset->n; t++) {
	    if (na(dataset->Z[1][t])) {
		miss = 1;
	    } else {
		empty = 0;
	    }
	}
    }

    if (empty) {
	infobox_printf(_("Warning: series %s is empty"), dataset->varname[1]);
    } else if (miss) {
	infobox(_("Warning: there were missing observations"));
    }

    register_data(DATA_APPENDED);
}

static gint maybe_exit_sheet (GtkWidget *w, Spreadsheet *sheet)
{
    if (sheet_is_modified(sheet)) {
	int resp = yes_no_dialog("gretl", _("Save changes?"),
				 sheet->win);

	if (resp == GRETL_YES) {
	    get_data_from_sheet(NULL, sheet);
	} else if (resp == GRETL_CANCEL) {
	    return FALSE;
	}
    }

    if (sheet->cmd == SHEET_NEW_DATASET) {
	empty_dataset_guard();
    }

    gtk_widget_destroy(sheet->win);

    return FALSE;
}

static gint sheet_delete_event (GtkWidget *w, GdkEvent *event,
				Spreadsheet *sheet)
{
    int resp;

    if (sheet_is_modified(sheet)) {
	resp = yes_no_cancel_dialog("gretl",
				    (sheet->matrix != NULL)? _("Save changes?") :
				    _("Do you want to save changes you have\n"
				      "made to the current data set?"),
				    sheet->win);
	if (resp == GRETL_YES) {
	    get_data_from_sheet(NULL, sheet);
	} else if (resp == GRETL_CANCEL) {
	    return TRUE;
	}
    }

    if (sheet->cmd == SHEET_NEW_DATASET) {
	empty_dataset_guard();
    }

    return FALSE;
}

static gint simple_exit_sheet (GtkWidget *w, Spreadsheet *sheet)
{
    gtk_widget_destroy(sheet->win);
    return FALSE;
}

static void size_matrix_window (Spreadsheet *sheet)
{
    int nc = sheet->datacols;
    int w, h;

    if (nc < 2) {
	nc = 2;
    }

    w = get_row_label_width(sheet) + nc * get_data_col_width() + 30;
    if (w > 640) {
	w = 640;
    }

    h = sheet->datarows * 20 + 160;
    if (h > 480) {
	h = 480;
    }

    gtk_window_resize(GTK_WINDOW(sheet->win), w, h);
}

static void size_scalars_window (Spreadsheet *sheet)
{
    int nc = 2;
    int w, h;

    w = get_row_label_width(sheet) + nc * get_data_col_width() +
	get_delete_col_width() + 30;
    if (w > 640) {
	w = 640;
    }

    h = sheet->datarows * 20 + 160;
    if (h > 480) {
	h = 480;
    }

    gtk_window_resize(GTK_WINDOW(sheet->win), w, h);
}

static void size_data_window (Spreadsheet *sheet, int hscroll)
{
    int w, h = 400;

    if (hscroll) {
	int nc = (sheet->varlist[0] > 4)? 4 : sheet->varlist[0];
	int ocw = get_obs_col_width();
	int dcw = get_data_col_width();
	int extra = 40;

	if (nc < 2) {
	    extra += 40;
	}
	w = ocw + nc * dcw + extra;
    } else {
	w = -1;
    }

    gtk_window_set_default_size(GTK_WINDOW(sheet->win), w, h);
}

/* hack to avoid losing a not-yet-committed edit to a cell
   in the sheet, on choosing Save, Apply, OK, etc */

static gboolean
button_entered (GtkWidget *w, GdkEventCrossing *e, Spreadsheet *sheet)
{
    if (sheet->entry != NULL) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(sheet->entry));
	GtkTreePath *path;
	gchar *pathstr;

	gtk_tree_view_get_cursor(GTK_TREE_VIEW(sheet->view), &path, NULL);

	if (path != NULL) {
	    pathstr = gtk_tree_path_to_string(path);
	    sheet->next = NEXT_SAME;
	    sheet_cell_edited(NULL, pathstr, s, sheet);
	    g_free(pathstr);
	    gtk_tree_path_free(path);
	}
    }

    return FALSE;
}

static void adjust_add_menu_state (Spreadsheet *sheet)
{
    sheet->flags |= (SHEET_ADD_OBS_OK | SHEET_INSERT_OBS_OK);

    if (complex_subsampled() || dataset->t2 < dataset->n - 1) {
	sheet->flags &= ~SHEET_ADD_OBS_OK;
	sheet->flags &= ~SHEET_INSERT_OBS_OK;
    } else if ((sheet->flags & SHEET_SHORT_VARLIST) ||
	       dataset_is_panel(dataset)) {
	sheet->flags &= ~SHEET_INSERT_OBS_OK;
    }
}

static void series_sheet_add_locator (Spreadsheet *sheet,
				      GtkWidget *hbox)
{
    GtkWidget *vbox = gtk_vbox_new(FALSE, 1);
    GtkWidget *status_box = gtk_hbox_new(FALSE, 1);
    gint w = get_obs_col_width();

    gtk_container_set_border_width(GTK_CONTAINER(status_box), 0);
    gtk_box_pack_start(GTK_BOX(vbox), status_box, FALSE, FALSE, 0);

    sheet->locator = gtk_statusbar_new();
    gtk_widget_set_size_request(sheet->locator, 2 * w, 20);
#if GTK_MAJOR_VERSION < 3
    gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(sheet->locator), FALSE);
#endif
    sheet->cid = gtk_statusbar_get_context_id(GTK_STATUSBAR(sheet->locator),
					      "current row and column");
    gtk_box_pack_start(GTK_BOX(status_box), sheet->locator, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);
}

static void sheet_add_toolbar (Spreadsheet *sheet, GtkWidget *vbox)
{
    GtkWidget *hbox, *tbar;
    GretlToolItem *items;
    GretlToolItem *item;
    GtkWidget *button;
    int i, n_items;

    hbox = gtk_hbox_new(FALSE, 0);
    tbar = gretl_toolbar_new(NULL);

    if (editing_scalars(sheet)) {
	items = scalar_items;
	n_items = n_scalar_items;
    } else {
	items = series_items;
	n_items = n_series_items;
    }

    for (i=0; i<n_items; i++) {
	item = &items[i];
	button = gretl_toolbar_insert(tbar, item, item->func, sheet, -1);
	if (!editing_scalars(sheet)) {
	    if (item->flag == SHEET_APPLY_BTN) {
		sheet->apply = button;
		gtk_widget_set_sensitive(sheet->apply, FALSE);
	    }
	    if (item->flag != SHEET_ADD_BTN) {
		g_signal_connect(G_OBJECT(button), "enter-notify-event",
				 G_CALLBACK(button_entered), sheet);
	    }
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), tbar, FALSE, FALSE, 0);

    if (!editing_scalars(sheet)) {
	series_sheet_add_locator(sheet, hbox);
    }

    window_add_winlist(sheet->win, hbox);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    if (!editing_scalars(sheet)) {
	adjust_add_menu_state(sheet);
    }
}

static void sheet_add_matrix_menu (Spreadsheet *sheet, GtkWidget *vbox)
{
    GtkActionGroup *actions;
    GtkWidget *mbar;

    sheet->ui = gtk_ui_manager_new();
    actions = gtk_action_group_new("SheetActions");
    gtk_action_group_set_translation_domain(actions, "gretl");

    gtk_action_group_add_actions(actions, matrix_items,
				 sizeof matrix_items / sizeof matrix_items[0],
				 sheet);
    gtk_ui_manager_add_ui_from_string(sheet->ui, matrix_ui, -1, NULL);

    gtk_ui_manager_insert_action_group(sheet->ui, actions, 0);
    g_object_unref(actions);

    mbar = gtk_ui_manager_get_widget(sheet->ui, "/menubar");
    gtk_box_pack_start(GTK_BOX(vbox), mbar, FALSE, FALSE, 0);
}

static void sheet_add_matrix_locator (Spreadsheet *sheet, GtkWidget *vbox)
{
    GtkWidget *status_box = gtk_hbox_new(FALSE, 1);
    gint w = get_row_label_width(NULL);

    gtk_container_set_border_width(GTK_CONTAINER(status_box), 0);
    gtk_box_pack_start(GTK_BOX(vbox), status_box, FALSE, FALSE, 0);

    sheet->locator = gtk_statusbar_new();
    gtk_widget_set_size_request(sheet->locator, 2 * w, 20);
#if GTK_MAJOR_VERSION < 3
    gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(sheet->locator), FALSE);
#endif
    sheet->cid = gtk_statusbar_get_context_id(GTK_STATUSBAR(sheet->locator),
					      "current row and column");
    gtk_box_pack_start(GTK_BOX(status_box), sheet->locator, FALSE, FALSE, 0);
}

/* this has to block when the user is defining a matrix in the
   course of responding to the function call dialog; otherwise
   it should not block
*/

static void real_show_spreadsheet (Spreadsheet **psheet, SheetCmd c,
				   int block)
{
    Spreadsheet *sheet = *psheet;
    GtkWidget *tmp, *scroller, *main_vbox;
    int hscroll = 1;
    int err = 0;

    sheet->win = gretl_gtk_window();

    if (sheet->matrix != NULL) {
	if (sheet->mname[0] != '\0') {
	    gchar *tmp = g_strdup_printf("gretl: %s", sheet->mname);

	    gtk_window_set_title(GTK_WINDOW(sheet->win), tmp);
	    g_free(tmp);
	} else {
	    gtk_window_set_title(GTK_WINDOW(sheet->win), _("gretl: edit matrix"));
	}
    } else if (c == SHEET_EDIT_SCALARS) {
	gtk_window_set_title(GTK_WINDOW(sheet->win), _("gretl: scalars"));
    } else {
	gtk_window_set_title(GTK_WINDOW(sheet->win), _("gretl: edit data"));
    }

    if (sheet->matrix != NULL) {
	size_matrix_window(sheet);
    } else if (c == SHEET_EDIT_SCALARS) {
	size_scalars_window(sheet);
	hscroll = 0;
    } else {
	if (sheet->varlist != NULL && sheet->varlist[0] < 3) {
	    /* make columns full "natural" width */
	    hscroll = 0;
	}
	size_data_window(sheet, hscroll);
    }

    if (block) {
	g_signal_connect(G_OBJECT(sheet->win), "destroy",
			 G_CALLBACK(gtk_main_quit), NULL);
    }

    g_signal_connect(G_OBJECT(sheet->win), "delete-event",
		     G_CALLBACK(sheet_delete_event), sheet);

    main_vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 5);
    gtk_container_add(GTK_CONTAINER(sheet->win), main_vbox);

    if (sheet->matrix != NULL) {
	sheet_add_matrix_menu(sheet, main_vbox);
	sheet_add_matrix_locator(sheet, main_vbox);
    } else {
	if (c != SHEET_EDIT_SCALARS) {
	    if (sheet->varlist != NULL) {
		int i;

		sheet->datacols = 0;

		for (i=sheet->varlist[0]; i>0; i--) {
		    if (series_is_hidden(dataset, sheet->varlist[i])) {
			gretl_list_delete_at_pos(sheet->varlist, i);
		    } else {
			sheet->datacols += 1;
		    }
		}
		if (sheet->varlist[0] < dataset->v - 1) {
		    sheet->flags |= SHEET_SHORT_VARLIST;
		}
	    }
	}
	sheet_add_toolbar(sheet, main_vbox);
    }

    gtk_widget_show_all(main_vbox);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   hscroll ? GTK_POLICY_AUTOMATIC :
				   GTK_POLICY_NEVER,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroller),
					GTK_SHADOW_IN);

    build_sheet_view(sheet);
    gtk_container_add(GTK_CONTAINER(scroller), sheet->view);
    gtk_box_pack_start(GTK_BOX(main_vbox), scroller, TRUE, TRUE, TRUE);

    gtk_widget_show(sheet->view);
    gtk_widget_show(scroller);

    if (sheet->matrix != NULL) {
	g_signal_connect(G_OBJECT(sheet->win), "destroy",
			 G_CALLBACK(free_matrix_sheet), sheet);
    } else if (c == SHEET_EDIT_SCALARS) {
	g_signal_connect(G_OBJECT(sheet->win), "destroy",
			 G_CALLBACK(free_spreadsheet), psheet); /* FIXME? */
    } else {
	g_signal_connect(G_OBJECT(sheet->win), "destroy",
			 G_CALLBACK(free_spreadsheet), psheet);
    }

    if (sheet->matrix != NULL) {
	/* control buttons (for edting matrices only) */
	GtkWidget *button_box;

	button_box = gtk_hbutton_box_new();
	gtk_button_box_set_layout(GTK_BUTTON_BOX(button_box),
				  GTK_BUTTONBOX_END);
	gtk_box_set_spacing(GTK_BOX(button_box), 10);
	gtk_widget_show(button_box);
	gtk_box_pack_start(GTK_BOX(main_vbox), button_box, FALSE, FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(button_box), 0);

	if (sheet->matrix != NULL && sheet->oldmat == NULL) {
	    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
	    gtk_container_add(GTK_CONTAINER(button_box), tmp);
	    g_signal_connect(G_OBJECT(tmp), "enter-notify-event",
			     G_CALLBACK(button_entered), sheet);
	    g_signal_connect(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(get_data_from_sheet), sheet);
	    g_signal_connect(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(simple_exit_sheet), sheet);
	    gtk_widget_show(tmp);
	} else if (sheet->matrix != NULL) {
	    /* editing a new matrix */
	    tmp = gtk_button_new_from_stock(GTK_STOCK_SAVE_AS);
	    gtk_container_add(GTK_CONTAINER(button_box), tmp);
	    g_signal_connect(G_OBJECT(tmp), "enter-notify-event",
			     G_CALLBACK(button_entered), sheet);
	    g_signal_connect(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(matrix_save_as), sheet);
	    gtk_widget_show(tmp);

	    tmp = gtk_button_new_from_stock(GTK_STOCK_SAVE);
	    gtk_container_add(GTK_CONTAINER(button_box), tmp);
	    g_signal_connect(G_OBJECT(tmp), "enter-notify-event",
			     G_CALLBACK(button_entered), sheet);
	    g_signal_connect(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(get_data_from_sheet), sheet);
	    gtk_widget_show(tmp);
	    sheet->save = tmp;
	    gtk_widget_set_sensitive(sheet->save, FALSE);

	    tmp = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
	    gtk_container_add(GTK_CONTAINER(button_box), tmp);
	    g_signal_connect(G_OBJECT(tmp), "enter-notify-event",
			     G_CALLBACK(button_entered), sheet);
	    g_signal_connect(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(maybe_exit_sheet), sheet);
	    gtk_widget_show(tmp);
	} else {
	    /* doing scalars */
	    tmp = gtk_button_new_from_stock(GTK_STOCK_ADD);
	    gtk_container_add(GTK_CONTAINER(button_box), tmp);
	    g_signal_connect(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(add_scalar_callback), sheet);
	    gtk_widget_show(tmp);

	    tmp = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
	    gtk_container_add(GTK_CONTAINER(button_box), tmp);
	    g_signal_connect(G_OBJECT(tmp), "enter-notify-event",
			     G_CALLBACK(button_entered), sheet);
	    g_signal_connect(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(maybe_exit_sheet), sheet);
	    gtk_widget_show(tmp);

	    set_scalar_edit_callback(scalars_changed_callback);
	}
    }

    if (c == SHEET_EDIT_MATRIX) {
	err = add_matrix_data_to_sheet(sheet);
    } else if (c == SHEET_EDIT_SCALARS) {
	err = add_scalars_to_sheet(sheet);
    } else {
	err = add_data_to_sheet(sheet, c);
    }

    if (err) {
	gtk_widget_destroy(sheet->win);
	return;
    }

    g_signal_connect(G_OBJECT(sheet->view), "key-press-event",
		     G_CALLBACK(catch_spreadsheet_key), sheet);
    g_signal_connect(G_OBJECT(sheet->view), "button-press-event",
		     G_CALLBACK(catch_spreadsheet_click),
		     sheet);

    window_list_add(sheet->win, SSHEET);
    gtk_widget_show(sheet->win);
    select_first_editable_cell(sheet);
    gtk_widget_grab_focus(sheet->view);

    if (editing_series(sheet)) {
	/* we can't have the user making confounding changes elsewhere,
	   while editing the dataset here
	*/
	set_dataset_locked(TRUE);
    }

    if (block) {
	gretl_set_window_modal(sheet->win);
	gtk_main();
    }
}

void show_spreadsheet (SheetCmd c)
{
    static Spreadsheet *sheet;

    if (dataset->v == 1) {
	warnbox(_("Please add a variable to the dataset first"));
	return;
    }

    if (sheet != NULL) {
	gtk_window_present(GTK_WINDOW(sheet->win));
	return;
    }

    sheet = spreadsheet_new(c, 0);
    if (sheet == NULL) {
	return;
    }

    real_show_spreadsheet(&sheet, c, 0);
}

void show_spreadsheet_for_series (int varnum)
{
    static Spreadsheet *sheet;
    int c = SHEET_EDIT_VARLIST;

    if (sheet != NULL) {
	/* FIXME? */
	gtk_window_present(GTK_WINDOW(sheet->win));
	return;
    }

    sheet = spreadsheet_new(c, varnum);
    if (sheet == NULL) {
	return;
    }

    real_show_spreadsheet(&sheet, c, 0);
}

void edit_scalars (void)
{
    static Spreadsheet *sheet;

    if (sheet != NULL) {
	maybe_update_scalars_sheet(sheet);
	gtk_window_present(GTK_WINDOW(sheet->win));
	return;
    }

    sheet = spreadsheet_new(SHEET_EDIT_SCALARS, 0);
    if (sheet == NULL) {
	return;
    }

    scalars_sheet = sheet;

    sheet->datarows = n_user_scalars();
    sheet->datacols = 1;

    real_show_spreadsheet(&sheet, SHEET_EDIT_SCALARS, 0);
}

void sync_scalars_window (void)
{
    if (scalars_sheet != NULL) {
	maybe_update_scalars_sheet(scalars_sheet);
    }
}

struct gui_matrix_spec {
    char name[VNAMELEN];
    int rows;
    int cols;
    double fill;
    char *formula;
    int uselist;
    gretl_matrix *m;
};

struct mdialog {
    GtkWidget *dlg;
    GtkWidget *nentry;
    GtkWidget *ventry;
    GtkWidget *numerics;
    GtkWidget *formula;
    struct gui_matrix_spec *spec;
    int *retp;
};

static void spin_call (GtkWidget *s, int *n)
{
    *n = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(s));
}

static void matrix_dialog_ok (GtkWidget *w, struct mdialog *mdlg)
{
    const char *etxt;
    int err;

    etxt = gtk_entry_get_text(GTK_ENTRY(mdlg->nentry));

    if (etxt == NULL || *etxt == '\0') {
	infobox(_("You must give the matrix a name"));
	*mdlg->retp = GRETL_CANCEL;
	return;
    }

    if (gtk_widget_is_sensitive(mdlg->nentry)) {
	err = gui_validate_varname(etxt,
				   GRETL_TYPE_MATRIX,
				   mdlg->dlg);
	if (err) {
	    *mdlg->retp = GRETL_CANCEL;
	    return;
	}
    }

    strcpy(mdlg->spec->name, etxt);
    mdlg->spec->uselist = 0;

    if (mdlg->formula != NULL &&
	gtk_widget_is_sensitive(mdlg->formula)) {
	etxt = gtk_entry_get_text(GTK_ENTRY(mdlg->formula));
	if (etxt == NULL || *etxt == '\0') {
	    errbox(_("The matrix formula is empty"));
	    *mdlg->retp = GRETL_CANCEL;
	    return;
	}
	mdlg->spec->formula = g_strdup(etxt);
    } else if (mdlg->numerics != NULL &&
	       gtk_widget_is_sensitive(mdlg->numerics)) {
	double x;

	etxt = gtk_entry_get_text(GTK_ENTRY(mdlg->ventry));
	x = gui_double_from_string(etxt, &err);
	if (err) {
	    gtk_editable_select_region(GTK_EDITABLE(mdlg->ventry), 0, -1);
	    *mdlg->retp = GRETL_CANCEL;
	    return;
	}
	mdlg->spec->fill = x;
    } else {
	/* will select data series */
	mdlg->spec->uselist = 1;
    }

    gtk_widget_destroy(mdlg->dlg);
}

static gint choose_series (GtkWidget *w, struct mdialog *mdlg)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	if (mdlg->numerics != NULL) {
	    gtk_widget_set_sensitive(mdlg->numerics, FALSE);
	}
	if (mdlg->formula != NULL) {
	    gtk_widget_set_sensitive(mdlg->formula, FALSE);
	}
    }

    return FALSE;
}

static gint choose_formula (GtkWidget *w, struct mdialog *mdlg)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	if (mdlg->numerics != NULL) {
	    gtk_widget_set_sensitive(mdlg->numerics, FALSE);
	}
	if (mdlg->formula != NULL) {
	    gtk_widget_set_sensitive(mdlg->formula, TRUE);
	}
    }

    return FALSE;
}

static gint choose_numeric (GtkWidget *w, struct mdialog *mdlg)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	if (mdlg->numerics != NULL) {
	    gtk_widget_set_sensitive(mdlg->numerics, TRUE);
	}
	if (mdlg->formula != NULL) {
	    gtk_widget_set_sensitive(mdlg->formula, FALSE);
	}
    }

    return FALSE;
}

static int new_matrix_dialog (struct gui_matrix_spec *spec,
			      GtkWidget *parent,
			      int fncall)
{
    struct mdialog mdlg;
    GSList *group = NULL;
    GtkWidget *dlg;
    GtkWidget *hbox;
    GtkWidget *vbox;
    GtkWidget *tab;
    GtkWidget *rb;
    GtkWidget *w;
    gchar *pname = NULL;
    int maxdim = 1000;
    int series_ok = 1;
    int ret = GRETL_CANCEL;

    dlg = gretl_dialog_new("Matrix", parent,
			   GRETL_DLG_BLOCK | GRETL_DLG_MODAL);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    mdlg.dlg = dlg;
    mdlg.spec = spec;
    mdlg.numerics = NULL;
    mdlg.formula = NULL;
    mdlg.retp = &ret;

    if (fncall && parent != NULL) {
	/* we're being called from the function-call dialog
	   for a gfn package: pick up some relevant info
	*/
	get_fncall_param_info(parent, &series_ok, &pname);
    }

    /* top label */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("New matrix"));
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* matrix name entry */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Name:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    w = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(w), 24);
    if (*spec->name != '\0') {
	gtk_entry_set_text(GTK_ENTRY(w), spec->name);
	gtk_widget_set_sensitive(w, FALSE);
    } else {
	if (pname != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(w), pname);
	}
	gtk_entry_set_activates_default(GTK_ENTRY(w), TRUE);
    }
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    mdlg.nentry = w;

    g_free(pname);

    /* matrix construction options */

    if (series_ok) {
	/* option: build from series */
	hbox = gtk_hbox_new(FALSE, 5);
	rb = gtk_radio_button_new_with_label(NULL, _("Build from series"));
	g_signal_connect(G_OBJECT(rb), "clicked", G_CALLBACK(choose_series), &mdlg);
	gtk_box_pack_start(GTK_BOX(hbox), rb, TRUE, TRUE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rb));
    }

    /* option: build from formula */
    hbox = gtk_hbox_new(FALSE, 5);
    rb = gtk_radio_button_new_with_label(group, _("Build from formula"));
    g_signal_connect(G_OBJECT(rb), "clicked", G_CALLBACK(choose_formula), &mdlg);
    gtk_box_pack_start(GTK_BOX(hbox), rb, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(" ");
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    w = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(w), MAXLEN);
    gtk_entry_set_width_chars(GTK_ENTRY(w), 32);
    gtk_entry_set_activates_default(GTK_ENTRY(w), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    if (fncall) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rb), TRUE);
    } else {
	gtk_widget_set_sensitive(w, FALSE);
    }
    mdlg.formula = w;

    /* option: build numerically */
    hbox = gtk_hbox_new(FALSE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rb));
    rb = gtk_radio_button_new_with_label(group, _("Build numerically"));
    g_signal_connect(G_OBJECT(rb), "clicked", G_CALLBACK(choose_numeric), &mdlg);
    gtk_box_pack_start(GTK_BOX(hbox), rb, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tab = gtk_table_new(3, 2, FALSE);
    gtk_table_set_col_spacing(GTK_TABLE(tab), 0, 5);

    w = gtk_label_new(_("Number of rows:"));
    gtk_misc_set_alignment(GTK_MISC(w), 0, 1);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, 0, 1);
    w = gtk_spin_button_new_with_range(1, maxdim, 1);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(w), 0);
    gtk_entry_set_width_chars(GTK_ENTRY(w), 4);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), (gdouble) spec->rows);
    g_signal_connect(G_OBJECT(w), "value-changed",
		     G_CALLBACK(spin_call), &spec->rows);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, 0, 1);

    w = gtk_label_new(_("Number of columns:"));
    gtk_misc_set_alignment(GTK_MISC(w), 0, 1);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, 1, 2);
    w = gtk_spin_button_new_with_range(1, maxdim, 1);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(w), 0);
    gtk_entry_set_width_chars(GTK_ENTRY(w), 4);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), (gdouble) spec->cols);
    g_signal_connect(G_OBJECT(w), "value-changed",
		     G_CALLBACK(spin_call), &spec->cols);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, 1, 2);

    w = gtk_label_new(_("Initial fill value:"));
    gtk_misc_set_alignment(GTK_MISC(w), 0, 1);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, 2, 3);
    w = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(w), 16);
    gtk_entry_set_activates_default(GTK_ENTRY(w), TRUE);
    gtk_entry_set_text(GTK_ENTRY(w), "0");
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, 2, 3);
    mdlg.ventry = w;

    w = gtk_label_new(" ");
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), tab, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    if (series_ok || fncall) {
	gtk_widget_set_sensitive(hbox, FALSE);
    } else {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rb), TRUE);
    }

    mdlg.numerics = hbox;

    /* control buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    w = cancel_delete_button(hbox, dlg);
    w = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(w), "clicked", G_CALLBACK(matrix_dialog_ok), &mdlg);
    gtk_widget_grab_default(w);

    gtk_widget_show_all(dlg);

    return ret;
}

static int gui_matrix_from_list (selector *sr)
{
    struct gui_matrix_spec *s = selector_get_data(sr);
    const char *buf = selector_list(sr);
    int *list;
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	errbox("No variables are selected");
	return 1;
    } else {
	list = command_list_from_string(buf, &err);
    }

    if (err) {
	/* error message handled already */
	return err;
    }

    s->m = gretl_matrix_data_subset(list, dataset,
				    dataset->t1, dataset->t2,
				    M_MISSING_SKIP, &err);

    if (!err) {
	err = user_var_add_or_replace(s->name,
				      GRETL_TYPE_MATRIX,
				      s->m);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	/* record in command log */
	PRN *prn = NULL;

	if (bufopen(&prn) == 0) {
	    int i;

	    pprintf(prn, "matrix %s = {", s->name);
	    for (i=1; i<=list[0]; i++) {
		pprintf(prn, "%s", dataset->varname[list[i]]);
		if (i < list[0]) {
		    pputc(prn, ',');
		} else {
		    pputc(prn, '}');
		}
	    }
	    add_command_to_stack(gretl_print_get_buffer(prn), 0);
	    gretl_print_destroy(prn);
	}
    }

    free(list);

    return err;
}

static int matrix_from_formula (struct gui_matrix_spec *s)
{
    gchar *genline = g_strdup_printf("matrix %s = %s", s->name,
				     s->formula);
    int err;

    err = gui_run_genr(genline, dataset, OPT_NONE, NULL);

    if (err) {
	gui_errmsg(err);
    } else {
	s->m = get_matrix_by_name(s->name);
	if (s->m == NULL) {
	    err = 1;
	} else {
	    add_command_to_stack(genline, 0);
	}
    }

    g_free(genline);

    return err;
}

static int matrix_from_spec (struct gui_matrix_spec *s)
{
    gchar *genline;
    int err;

    if (na(s->fill)) {
	genline = g_strdup_printf("matrix %s = ones(%d,%d) * 0/0",
				  s->name, s->rows, s->cols);
    } else {
	genline = g_strdup_printf("matrix %s = zeros(%d,%d) + %.15g",
				  s->name, s->rows, s->cols, s->fill);
    }

    err = gui_run_genr(genline, dataset, OPT_NONE, NULL);

    if (err) {
	gui_errmsg(err);
    } else {
	s->m = get_matrix_by_name(s->name);
	if (s->m == NULL) {
	    err = 1;
	} else {
	    add_command_to_stack(genline, 0);
	}
    }

    g_free(genline);

    return err;
}

static void gui_matrix_spec_init (struct gui_matrix_spec *s,
				  gretl_matrix *m,
				  const char *name)
{
    *s->name = '\0';
    s->rows = 2;
    s->cols = 2;
    s->fill = 0;
    s->formula = NULL;
    s->uselist = 1;
    s->m = m;

    if (name != NULL) {
	strncat(s->name, name, VNAMELEN - 1);
    }
}

static void gui_edit_matrix (gretl_matrix *m, const char *name)
{
    Spreadsheet *sheet = NULL;
    int err = 0;

    if (m == NULL || name == NULL) {
	return;
    }

    sheet = spreadsheet_new(SHEET_EDIT_MATRIX, 0);
    if (sheet == NULL) {
	return;
    }

    strcpy(sheet->mname, name);
    sheet->oldmat = m;
    sheet->matrix = gretl_matrix_copy(m);
    if (sheet->matrix == NULL) {
	err = E_ALLOC;
    }

    if (err) {
	nomem();
	free(sheet);
	sheet = NULL;
    } else {
	sheet->colnames = gretl_matrix_get_colnames(m);
	sheet->rownames = gretl_matrix_get_rownames(m);
	sheet->datarows = gretl_matrix_rows(sheet->matrix);
	sheet->datacols = gretl_matrix_cols(sheet->matrix);
	real_show_spreadsheet(&sheet, SHEET_EDIT_MATRIX, 0);
    }

    if (sheet != NULL) {
	/* protect matrix from deletion while editing */
	user_var *u = get_user_var_by_name(name);

	if (u != NULL) {
	    g_object_set_data(G_OBJECT(sheet->win), "object", u);
	}
    }
}

static void edit_new_arg_matrix (gretl_matrix *m, const char *name)
{
    Spreadsheet *sheet = NULL;

    if (m == NULL || name == NULL) {
	return;
    }

    sheet = spreadsheet_new(SHEET_EDIT_MATRIX, 0);
    if (sheet == NULL) {
	return;
    }

    strcpy(sheet->mname, name);
    sheet->oldmat = m;
    sheet->matrix = gretl_matrix_copy(m);

    if (sheet->matrix == NULL) {
	nomem();
	free(sheet);
	sheet = NULL;
	return;
    }

    sheet->colnames = NULL;
    sheet->rownames = NULL;
    sheet->datarows = gretl_matrix_rows(sheet->matrix);
    sheet->datacols = gretl_matrix_cols(sheet->matrix);
    real_show_spreadsheet(&sheet, SHEET_EDIT_MATRIX, 1);
}

/* note that both @m and @name may be NULL depending on how
   we are called */

static void real_gui_new_matrix (gretl_matrix *m, const char *name,
				 GtkWidget *parent, int fncall)
{
    struct gui_matrix_spec spec;
    int resp;

    gui_matrix_spec_init(&spec, m, name);

    resp = new_matrix_dialog(&spec, parent, fncall);
    if (canceled(resp)) {
	return;
    }

    if (spec.uselist) {
	/* matrix from listed vars */
	simple_selection_with_data(DEFINE_MATRIX, _("Define matrix"),
				   gui_matrix_from_list,
				   NULL, &spec);
    } else if (spec.formula != NULL) {
	/* matrix from genr-style formula */
	matrix_from_formula(&spec);
	free(spec.formula);
    } else {
	/* numerical specification: open editor if OK */
	int err = matrix_from_spec(&spec);

	if (!err) {
	    if (fncall) {
		edit_new_arg_matrix(spec.m, spec.name);
	    } else {
		gui_edit_matrix(spec.m, spec.name);
	    }
	}
    }
}

/* callback for "Define matrix..." in main window, and
   also "+" to add matrix arg in fncall.c
*/

void gui_new_matrix (GtkWidget *parent)
{
    real_gui_new_matrix(NULL, NULL, parent, 0);
}

void fncall_add_matrix (GtkWidget *parent)
{
    real_gui_new_matrix(NULL, NULL, parent, 1);
}

static void display_complex_matrix (gretl_matrix *m,
				    const char *name)
{
    gchar *title = NULL;
    PRN *prn = NULL;
    int hsize, vsize;

    if (bufopen(&prn)) {
	return;
    }

    hsize = m->cols * 10;
    if (hsize < 70) {
	hsize = 70;
    } else if (hsize > 120) {
	hsize = 120;
    }

    vsize = 24 * (m->rows + 2);
    if (vsize < 200) {
	vsize = 200;
    } else if (vsize > 500) {
	vsize = 500;
    }

    gretl_cmatrix_print(m, name, prn);
    title = g_strdup_printf("gretl: %s %s", _("complex matrix"), name);
    preset_viewer_flag(VWIN_NO_SAVE);
    view_buffer(prn, hsize, vsize, title, PRINT, m);
    g_free(title);
}

static void display_real_submatrix (gretl_matrix *m,
				    const char *name)
{
    gretl_matrix *tmp;
    gchar *title = NULL;
    PRN *prn = NULL;
    double mij;
    int hsize, vsize;
    int show = 50;
    int r, c;
    int i, j;

    if (bufopen(&prn)) {
	return;
    }

    /* We come here only if the matrix is too big to put
       into an editable window (r * c > 36000).
    */

    if (m->rows > 50 && m->cols > 50) {
	r = show;
	c = show;
    } else if (m->rows == 1) {
	r = 1;
	c = show;
    } else if (m->cols == 1) {
	c = 1;
	r = show;
    } else if (m->rows > show) {
	r = show;
	c = MIN(show, m->cols);
    } else {
	c = show;
	r = MIN(show, m->rows);
    }

    hsize = c * 10;
    if (hsize < 70) {
	hsize = 70;
    } else if (hsize > 120) {
	hsize = 120;
    }

    vsize = 24 * (r + 2);
    if (vsize < 200) {
	vsize = 200;
    } else if (vsize > 500) {
	vsize = 500;
    }

    if (r == 1 || c == 1) {
	pprintf(prn, "%s is %d x %d; showing the first %d elements\n\n",
		name, m->rows, m->cols, r == 1 ? c : r);
    } else {
	pprintf(prn, "%s is %d x %d; showing north-west %d x %d submatrix\n\n",
		name, m->rows, m->cols, r, c);
    }
    tmp = gretl_matrix_alloc(r, c);
    for (j=0; j<c; j++) {
	for (i=0; i<r; i++) {
	    mij = gretl_matrix_get(m, i, j);
	    gretl_matrix_set(tmp, i, j, mij);
	}
    }

    gretl_matrix_print_to_prn(tmp, NULL, prn);
    title = g_strdup_printf("gretl: %s %s", _("matrix"), name);
    preset_viewer_flag(VWIN_NO_SAVE);
    view_buffer(prn, hsize, vsize, title, PRINT, NULL);
    g_free(title);
    gretl_matrix_free(tmp);
}

void edit_or_view_matrix (const char *name, GtkWidget *parent)
{
    user_var *u = get_user_var_by_name(name);
    gretl_matrix *m = user_var_get_value(u);
    GtkWidget *w;

    if (m == NULL) {
	errbox_printf(_("Couldn't open '%s'"), name);
	return;
    }

    /* do we already have a window open for this matrix? */

    if (m->is_complex) {
	windata_t *vwin = get_viewer_for_data(m);

	if (vwin != NULL) {
	    gtk_window_present(GTK_WINDOW(vwin->main));
	    return;
	}
    } else {
	w = get_window_for_data(u);
	if (w != NULL) {
	    gtk_window_present(GTK_WINDOW(w));
	    return;
	}
    }

    if (m->is_complex) {
	/* we'll just display it for now */
	display_complex_matrix(m, name);
    } else if (gretl_is_null_matrix(m)) {
	real_gui_new_matrix(m, name, parent, 0);
    } else {
	int n = m->rows * m->cols;

	if (n > 36000) {
	    display_real_submatrix(m, name);
	} else {
	    gui_edit_matrix(m, name);
	}
    }
}

/* mechanism for locking dataset against changes while editing */

static int locked;

void set_dataset_locked (gboolean s)
{
    locked = s;
    flip(mdata->ui, "/menubar/Data", !s);
    flip(mdata->ui, "/menubar/Sample", !s);
}

int dataset_locked (void)
{
    if (locked) {
	errbox(_("The dataset cannot be modified at present"));
    }

    return locked;
}

/* mechanism for reformatting displayed matrix */

struct format_adjuster {
    int digits;
    char fmtchar;
    char oldfmt;
    gboolean custom;
    GtkWidget *spin;
    GtkWidget *combo;
    Spreadsheet *sheet;
};

static void sheet_toggle_custom (GtkWidget *w, struct format_adjuster *fa)
{
    gboolean std = button_is_active(w);

    fa->custom = !std;
    gtk_widget_set_sensitive(fa->spin, fa->custom);
    gtk_widget_set_sensitive(fa->combo, fa->custom);
}

static void reformat_sheet_callback (GtkButton *b, struct format_adjuster *fa)
{
    Spreadsheet *sheet = fa->sheet;

    if (!(sheet->flags & SHEET_CUSTOM_FMT) && !fa->custom) {
	/* no change */
	return;
    }

    if (fa->custom) {
	fa->digits = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(fa->spin));
	fa->fmtchar = gtk_combo_box_get_active(GTK_COMBO_BOX(fa->combo)) == 0 ?
	    'g' : 'f';
	if (fa->digits == sheet->digits && fa->fmtchar == fa->oldfmt) {
	    /* no change */
	    return;
	} else {
	    sheet->flags |= SHEET_CUSTOM_FMT;
	    sheet->digits = fa->digits;
	    strcpy(sheet->numfmt, fa->fmtchar == 'g' ? "%.*g" : "%.*f");
	}
    } else {
	sheet->flags &= ~SHEET_CUSTOM_FMT;
	sheet->digits = DBL_DIG;
	strcpy(sheet->numfmt, "%.*g");
    }

    reformat_sheet_matrix(sheet);
}

static char lastchar (const char *s)
{
    int n = strlen(s);

    return s[n-1];
}

static void sheet_number_format_dialog (Spreadsheet *sheet)
{
    struct format_adjuster fa;
    GtkWidget *dlg, *vbox;
    GtkWidget *tmp, *hbox;
    GtkWidget *b1, *b2;
    GSList *group;
    int std = 1;

    dlg = gretl_dialog_new(_("gretl: data format"), sheet->win,
			   GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    if (sheet->flags & SHEET_CUSTOM_FMT) {
	std = 0;
    }

    fa.digits = sheet->digits;
    fa.fmtchar = fa.oldfmt = lastchar(sheet->numfmt);
    fa.sheet = sheet;

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Select data format"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    /* standard format radio option */
    hbox = gtk_hbox_new(FALSE, 5);
    b1 = gtk_radio_button_new_with_label(NULL, _("Standard format"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), std);
    gtk_box_pack_start(GTK_BOX(hbox), b1, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    /* custom format radio option */
    hbox = gtk_hbox_new(FALSE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Show"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), !std);
    gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 5);

    /* with spin button for number of digits */
    fa.spin = gtk_spin_button_new_with_range(1, 16, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(fa.spin), sheet->digits);
    gtk_widget_set_sensitive(fa.spin, !std);
    gtk_box_pack_start(GTK_BOX(hbox), fa.spin, FALSE, FALSE, 0);

    /* and selector for digits / decimal places */
    fa.combo = gtk_combo_box_text_new();
    combo_box_append_text(fa.combo, _("significant figures"));
    combo_box_append_text(fa.combo, _("decimal places"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(fa.combo),
			     fa.fmtchar == 'g' ? 0 : 1);
    gtk_widget_set_sensitive(fa.combo, !std);
    gtk_box_pack_start(GTK_BOX(hbox), fa.combo, FALSE, FALSE, 5);

    /* connect toggle signal */
    g_signal_connect(G_OBJECT(b1), "toggled",
                     G_CALLBACK(sheet_toggle_custom), &fa);

    /* pack the custom line */
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    cancel_delete_button(hbox, dlg);

    /* OK button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(reformat_sheet_callback), &fa);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dlg);
}

static void matrix_format_callback (GtkAction *action, gpointer data)
{
    Spreadsheet *sheet = (Spreadsheet *) data;

    sheet_number_format_dialog(sheet);
}
