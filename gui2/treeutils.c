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

/* treeutils.c for gretl */

#include "gretl.h"
#include "treeutils.h"
#include "datafiles.h"
#include "dlgutils.h"

/* special comparator which always preserves "const" (variable 0) in
   the first position when sorting variables by name */

static gint list_alpha_compare (GtkTreeModel *model, 
				GtkTreeIter *a, GtkTreeIter *b,
				gpointer p)
{
    gchar *str_a, *str_b;
    gchar *id_a, *id_b;
    GtkSortType order;
    gint scol, ret;

    gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model),
					 &scol, &order);

    gtk_tree_model_get(model, a, 0, &id_a, scol, &str_a, -1);
    gtk_tree_model_get(model, b, 0, &id_b, scol, &str_b, -1);

    if (atoi(id_a) == 0) {
	ret = (order == GTK_SORT_DESCENDING)? 1 : 0;
    } else if (atoi(id_b) == 0) {
	ret = (order == GTK_SORT_DESCENDING)? 0 : 1;
    } else {
	ret = strcmp(str_a, str_b);
    }    

    g_free(id_a);
    g_free(id_b);

    g_free(str_a);
    g_free(str_b);
    
    return ret;
}

/* special comparator which preserves 0 in first position when sorting
   variables by ID number */

static gint list_id_compare (GtkTreeModel *model, 
			     GtkTreeIter *a, GtkTreeIter *b,
			     gpointer p)
{
    gchar *id_a, *id_b;
    GtkSortType order;
    gint ia, ib, scol, ret;

    gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model),
					 &scol, &order);
    
    gtk_tree_model_get(model, a, 0, &id_a, -1);
    gtk_tree_model_get(model, b, 0, &id_b, -1);

    ia = atoi(id_a);
    ib = atoi(id_b);

    if (ia == 0) {
	ret = (order == GTK_SORT_DESCENDING)? 1 : 0;
    } else if (ib == 0) {
	ret = (order == GTK_SORT_DESCENDING)? 0 : 1;
    } else {
	ret = ia - ib;
    }    

    g_free(id_a);
    g_free(id_b);
    
    return ret;
}

/* special comparator to put the series in a database window back into
   "natural" order (by position in the database) in place of a
   descending sort action in the first (series name) column
*/

static gint db_series_compare (GtkTreeModel *model, 
			       GtkTreeIter *a, GtkTreeIter *b,
			       GtkTreeViewColumn *col)
{
    GtkSortType order;
    gint scol, ret;

    gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model),
					 &scol, &order);

    if (order == GTK_SORT_ASCENDING) {
	/* alphabetize as usual */
	gchar *sa, *sb;

	gtk_tree_model_get(model, a, 0, &sa, -1);
	gtk_tree_model_get(model, b, 0, &sb, -1);
	ret = strcmp(sa, sb);
	g_free(sa);
	g_free(sb);
    } else {
	/* 'natural' order in place of descending sort:
	   use the last (integer) column */
	gint ia, ib;

	gtk_tree_model_get(model, a, 3, &ia, -1);
	gtk_tree_model_get(model, b, 3, &ib, -1);
	ret = ib - ia;
	gtk_tree_view_column_set_sort_indicator(col, FALSE);
    }
    
    return ret;
}

static gint 
compare_treenames (GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b,
		   gpointer p)
{
    gchar *s1, *s2;
    gint ret;

    gtk_tree_model_get(model, a, 0, &s1, -1);
    gtk_tree_model_get(model, b, 0, &s2, -1);

    ret = strcmp(gretl_lower(s1), gretl_lower(s2));

    g_free(s1);
    g_free(s2);

    return ret;    
}

/* sort in ascending order by the strings in the first column
   of the treeview in @vwin
*/

void presort_treelist (windata_t *vwin)
{
    GtkTreeModel *model;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));

    gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(model), 
					 0, GTK_SORT_ASCENDING);
    gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(model), 0,
				    compare_treenames, NULL, NULL);
}


static gboolean no_select_row_zero (GtkTreeSelection *selection,
				    GtkTreeModel *model,
				    GtkTreePath *path,
				    gboolean path_currently_selected,
				    gpointer data)
{
    return tree_path_get_row_number(path) != 0;
}

static void get_selected_varnum (GtkTreeModel *model, GtkTreePath *path,
				 GtkTreeIter *iter, int *v)
{
    gchar *id;

    gtk_tree_model_get(model, iter, 0, &id, -1);  
    *v = atoi(id);
    g_free(id);
}

int tree_selection_count (GtkTreeSelection *select, int *vnum)
{
    int selcount = 0;

    if (select != NULL) {
	selcount = gtk_tree_selection_count_selected_rows(select);
    }
    
    if (vnum != NULL && selcount == 1) {
	gtk_tree_selection_selected_foreach(select, 
					    (GtkTreeSelectionForeachFunc) 
					    get_selected_varnum,
					    vnum);	
    }

    return selcount;
}

int vwin_selection_count (windata_t *vwin, int *row)
{
    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
    GtkTreeSelection *sel = gtk_tree_view_get_selection(view);

    return tree_selection_count(sel, row);
}

static void my_gtk_entry_append_text (GtkEntry *entry, gchar *add)
{
    const gchar *old = gtk_entry_get_text(entry);
    gchar *new = g_strdup_printf("%s%s", old, add);
    
    gtk_entry_set_text(entry, new);
    gtk_editable_set_position(GTK_EDITABLE(entry), -1);
    g_free(new);
}

static void update_dialogs_from_varclick (int active_var)
{
    GtkWidget *active_edit_id = get_active_edit_id();
    GtkWidget *active_edit_name = get_active_edit_name();
    GtkWidget *active_edit_text = get_active_edit_text();
    const gchar *edttext;

    if (active_edit_id != NULL) {
	gchar addvar[9];

	edttext = gtk_entry_get_text(GTK_ENTRY(active_edit_id));
	if (*edttext != '\0') {
	    sprintf(addvar, " %d", active_var);
	} else {
	    sprintf(addvar, "%d", active_var);
	}
	my_gtk_entry_append_text(GTK_ENTRY(active_edit_id), addvar);
    } else if (active_edit_name != NULL) {
	edttext = gtk_entry_get_text(GTK_ENTRY(active_edit_name));
	my_gtk_entry_append_text(GTK_ENTRY(active_edit_name), 
				 dataset->varname[active_var]);
    } else if (active_edit_text != NULL) {
	GtkTextBuffer *tbuf;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(active_edit_text));
	if (tbuf != NULL) {
	    gtk_text_buffer_insert_at_cursor(tbuf, 
					     dataset->varname[active_var], 
					     -1);
	}
    }  
}

#ifdef OS_OSX

/* Support for command-click as replacement for Ctrl-click
   on OS X */

static gboolean maybe_do_meta_click (GdkEventButton *event, 
				     GtkTreeView *view,
				     GtkTreePath *path)
{
    if (event->state & GDK_MOD2_MASK) {
	/* The idea here is that we add to the current selection
	   and return TRUE to block further processing of the
	   click, so as to implement selection of non-contiguous
	   rows in the GtkTreeView. It works on Linux using 
	   GDK_CONTROL_MASK as the criterion (i.e. faking the
	   standard behavior of Ctrl-click as test of concept).
	*/
	GtkTreeSelection *sel = gtk_tree_view_get_selection(view);

	gtk_tree_selection_select_path(sel, path);
	return TRUE;
    } else {
	return FALSE;
    }
}

#endif

/* Respond to a mouse click in gretl's main window. This callback
   is connected after the one that deals with right-clicking to
   handle popups, so here we can assume it's not a right-click.
*/

gboolean main_varclick (GtkWidget *widget, GdkEventButton *event,
			windata_t *vwin)
{
    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
    GtkTreePath *path;
    gboolean ret = FALSE;

    if (dataset == NULL || dataset->n == 0) {
	return TRUE;
    }

    if (gtk_tree_view_get_path_at_pos(view, event->x, event->y, &path, 
				      NULL, NULL, NULL)) {
	gint row = tree_path_get_row_number(path);

	if (row == 0) {
	    ret = TRUE;
	} else {
	    gchar *varnum;

	    g_object_set_data(G_OBJECT(vwin->listbox), "active_row",
			      GINT_TO_POINTER(row));
	    tree_view_get_string(view, row, 0, &varnum);
	    vwin->active_var = atoi(varnum);
	    g_free(varnum);
	    update_dialogs_from_varclick(vwin->active_var);
#ifdef OS_OSX
	    ret = maybe_do_meta_click(event, view, path);
#endif
	} 

	gtk_tree_path_free(path);
    } else {
	/* clicked below the lines representing variables */
	return FALSE;
    }

    return ret;
}

static void
bool_col_toggled (GtkCellRendererToggle *cell, gchar *path_str, windata_t *vwin)
{
    GtkTreeView *treeview = GTK_TREE_VIEW(vwin->listbox);
    GtkTreeModel *model = gtk_tree_view_get_model(treeview);
    GtkTreePath *path = gtk_tree_path_new_from_string(path_str);
    GtkTreeIter iter;
    gboolean val;
    gint col;

    col = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(vwin->main), "boolcol"));
    gtk_tree_model_get_iter(model, &iter, path);
    gtk_tree_model_get(model, &iter, col, &val, -1);

    if (val) {
	gtk_tree_path_free(path);
	return;
    }

    if (vwin->role == FUNC_FILES) {
	vwin->active_var = atoi(path_str);
	browser_load_func(NULL, vwin);
    } 

    gtk_list_store_set(GTK_LIST_STORE(model), &iter, col, TRUE, -1);
    gtk_tree_path_free(path);
}

static gint catch_listbox_key (GtkWidget *w, GdkEventKey *event, 
			       windata_t *vwin)
{
    int key = event->keyval;

#ifdef MAC_NATIVE
    if (key == GDK_w && cmd_key(event)) {
	if (vwin != mdata) {
	    gtk_widget_destroy(vwin->main);
	}
	return TRUE;
    }
#endif	

    if (key == GDK_q) { 
	/* Q = quit */
	if (vwin != mdata) {
	    gtk_widget_destroy(vwin->main);
	}
	return TRUE;
    } else if (key == GDK_f) {
	/* F = find */
	if (vwin == mdata && !data_status) {
	    return TRUE;
	}
	if (event->state & GDK_CONTROL_MASK) {
	    listbox_find(NULL, vwin);
	    return TRUE;
	}	
    } 

    return FALSE;
}

static void check_series_pd (GtkTreeModel *model, GtkTreePath *path,
			     GtkTreeIter *iter, windata_t *vwin)
{
    static int pd0;
    gchar *tmp = NULL;

    if (model == NULL) {
	/* reset */
	pd0 = 0;
	return;
    }

    gtk_tree_model_get(model, iter, 2, &tmp, -1);
    if (tmp != NULL) {
	if (pd0 == 0) {
	    pd0 = *tmp;
	} else if (*tmp != pd0) {
	    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
	    GtkTreeSelection *sel;

	    sel = gtk_tree_view_get_selection(view);
	    gtk_tree_selection_unselect_iter(sel, iter);
	}
	g_free(tmp);
    }
}

static void set_active_row (GtkTreeModel *model, GtkTreePath *path,
			    GtkTreeIter *iter, windata_t *vwin)
{
    vwin->active_var = tree_path_get_row_number(path);
}

static void check_db_series_selection (GtkTreeSelection *sel, 
				       windata_t *vwin)
{
    int nsel = gtk_tree_selection_count_selected_rows(sel);

    if (nsel == 1) {
	gtk_tree_selection_selected_foreach(sel, 
					    (GtkTreeSelectionForeachFunc) 
					    set_active_row, vwin);
    } else {
	check_series_pd(NULL, NULL, NULL, NULL);
	gtk_tree_selection_selected_foreach(sel, 
					    (GtkTreeSelectionForeachFunc) 
					    check_series_pd, vwin);
    }  
}

static void id_col_clicked (GtkTreeViewColumn *column, GtkWidget *view)
{
    GtkTreeModel *model;
    GtkSortType order;
    gint scol;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(view));
    gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model),
					 &scol, &order);
    if (order == GTK_SORT_ASCENDING) {
	gtk_tree_view_column_set_sort_indicator(column, FALSE);
    }
}

#define db_series_window(v) (v->role == NATIVE_SERIES || \
                             v->role == RATS_SERIES || \
                             v->role == PCGIVE_SERIES || \
                             v->role == REMOTE_SERIES)

void vwin_add_list_box (windata_t *vwin, GtkBox *box, 
			int ncols, gboolean hidden_col,
			GType *types, const char **titles,
			int tree) 
{
    GtkListStore *lstore = NULL;
    GtkTreeStore *tstore = NULL;
    GtkWidget *view, *scroller;
    GtkCellRenderer *renderer;
    GtkCellRenderer *bool_renderer = NULL;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    int i, viscols = ncols;

    if (hidden_col) {
	viscols--;
    }

    if (tree) {
	tstore = gtk_tree_store_newv(ncols, types);
	vwin->listbox = view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(tstore));
    } else {
	lstore = gtk_list_store_newv(ncols, types);
	vwin->listbox = view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(lstore));
    }

    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), TRUE);

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 0, NULL);

    for (i=0; i<viscols; i++) {
	if (types[i] == G_TYPE_BOOLEAN) {
	    bool_renderer = gtk_cell_renderer_toggle_new();
	    g_object_set_data(G_OBJECT(vwin->main), "boolcol", GINT_TO_POINTER(i));
	    g_signal_connect(G_OBJECT(bool_renderer), "toggled",
			     G_CALLBACK(bool_col_toggled), vwin);
	    column = gtk_tree_view_column_new_with_attributes(_(titles[i]),
							      bool_renderer,
							      "active", i,
							      NULL);
	    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column),
					    GTK_TREE_VIEW_COLUMN_FIXED);
	    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
	    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	} else {
	    column = gtk_tree_view_column_new_with_attributes(_(titles[i]),
							      renderer,
							      "text", i, 
							      NULL);
	    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);

	    if ((vwin->role == FUNC_FILES || vwin->role == REMOTE_FUNC_FILES) && i == 1) {
		; /* sorting functions by version number is not meaningful */
	    } else {
		gtk_tree_view_column_set_sort_column_id(GTK_TREE_VIEW_COLUMN(column), i);
	    }

	    if (vwin != mdata) {
		g_object_set(G_OBJECT(column), "resizable", TRUE, NULL);
	    }

	    /* special actions on first column */
	    if (i == 0) {
		if (vwin == mdata) {
		    g_signal_connect(G_OBJECT(column), "clicked",
				     G_CALLBACK(id_col_clicked), view);
		} else if (db_series_window(vwin)) {
		    gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(lstore), 0, 
						    (GtkTreeIterCompareFunc) 
						    db_series_compare,
						    column, NULL);
		}
	    }
	}	
    }

    if (hidden_col) {
	column = gtk_tree_view_column_new_with_attributes(NULL,
							  renderer,
							  "text", i, 
							  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	gtk_tree_view_column_set_visible(column, FALSE);
    }

    /* set the selection properties on the tree view */
    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));

    if (vwin == mdata) { 
	/* gretl main window */
	gtk_tree_selection_set_mode(select, GTK_SELECTION_MULTIPLE);
	gtk_tree_selection_set_select_function(select, 
					       (GtkTreeSelectionFunc)
					       no_select_row_zero,
					       NULL, NULL);
	gtk_widget_set_events(view, GDK_POINTER_MOTION_MASK 
			      | GDK_POINTER_MOTION_HINT_MASK);
        g_signal_connect(G_OBJECT(view), "motion-notify-event",
			 G_CALLBACK(listbox_drag), NULL);
    } else if (db_series_window(vwin)) {
	gtk_tree_selection_set_mode(select, GTK_SELECTION_MULTIPLE);
	g_signal_connect(G_OBJECT(select), "changed",
			 G_CALLBACK(check_db_series_selection),
			 vwin);
    } else {
	gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
	g_signal_connect(G_OBJECT(select), "changed",
			 G_CALLBACK(listbox_select_row),
			 vwin);
    }

    g_signal_connect(G_OBJECT(view), "key-press-event",
		     G_CALLBACK(catch_listbox_key), vwin);
    g_signal_connect(G_OBJECT(view), "button-press-event",
		     G_CALLBACK(listbox_double_click), vwin);

    /* set sort properties on the tree model */
    if (vwin == mdata && tstore != NULL) {
	gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(tstore), 0, 
					(GtkTreeIterCompareFunc) 
					list_id_compare,
					NULL, NULL);
	gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(tstore), 1, 
					(GtkTreeIterCompareFunc) 
					list_alpha_compare,
					NULL, NULL);
	gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(tstore), 2, 
					(GtkTreeIterCompareFunc) 
					list_alpha_compare,
					NULL, NULL);
    } 

    if (lstore != NULL) {
	g_object_unref(G_OBJECT(lstore));
    } else if (tstore != NULL) {
	g_object_unref(G_OBJECT(tstore));
    }

    scroller = gtk_scrolled_window_new(NULL, NULL);

    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroller),
					GTK_SHADOW_IN);

    gtk_container_add(GTK_CONTAINER(scroller), view);
    gtk_box_pack_start(box, scroller, TRUE, TRUE, TRUE);

    gtk_widget_show(view);
    gtk_widget_show(scroller);
}

void set_main_colheads_clickable (gboolean s)
{
    GtkTreeViewColumn *col;
    int i;

    for (i=0; i<3; i++) {
	col = gtk_tree_view_get_column(GTK_TREE_VIEW(mdata->listbox), i);
	gtk_tree_view_column_set_clickable(col, s);
    }

    if (!s) {
	gtk_widget_grab_focus(mdata->listbox);
    }
}

void tree_view_get_bool (GtkTreeView *view, int row, int col, gboolean *val)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *path;

    model = gtk_tree_view_get_model(view);
    path = g_strdup_printf("%d", row);
    if (gtk_tree_model_get_iter_from_string(model, &iter, path)) {
	gtk_tree_model_get(model, &iter, col, val, -1);
    }
    g_free(path);
}

void tree_view_get_string (GtkTreeView *view, int row, int col, gchar **val)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *path;

    model = gtk_tree_view_get_model(view);
    path = g_strdup_printf("%d", row);
    if (gtk_tree_model_get_iter_from_string(model, &iter, path)) {
	gtk_tree_model_get(model, &iter, col, val, -1);
    }
    g_free(path);
}

static void tree_view_set_string (GtkTreeView *view, int row, int col, 
				  const gchar *val, int tstore)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *path;

    model = gtk_tree_view_get_model(view);
    path = g_strdup_printf("%d", row);
    if (gtk_tree_model_get_iter_from_string(model, &iter, path)) {
	if (tstore == 1) {
	    gtk_tree_store_set(GTK_TREE_STORE(model), &iter, col, val, -1);
	} else {
	    gtk_list_store_set(GTK_LIST_STORE(model), &iter, col, val, -1);
	}
    }
    g_free(path);
}

void list_store_set_string (GtkTreeView *view, int row, int col, const gchar *val)
{
    tree_view_set_string(view, row, col, val, 0);
}

void tree_store_set_string (GtkTreeView *view, int row, int col, const gchar *val)
{
    tree_view_set_string(view, row, col, val, 1);
}

void tree_view_get_int (GtkTreeView *view, int row, int col, int *val)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *path;

    model = gtk_tree_view_get_model(view);
    path = g_strdup_printf("%d", row);
    if (gtk_tree_model_get_iter_from_string(model, &iter, path)) {
	gtk_tree_model_get(model, &iter, col, val, -1);
    }
    g_free(path);
}

int tree_path_get_row_number (GtkTreePath *path)
{
    return gtk_tree_path_get_indices(path)[0];
}

void tree_model_get_iter_last (GtkTreeModel *mod, GtkTreeIter *iter)
{
    GtkTreeIter first;
    gboolean s;

    s = gtk_tree_model_get_iter_first(mod, &first);
    *iter = first;
    while (s && gtk_tree_model_iter_next(mod, &first)) {
	*iter = first;
    }
}

gboolean tree_model_iter_prev (GtkTreeModel *mod, GtkTreeIter *iter)
{
    GtkTreePath *path;
    gboolean s;

    path = gtk_tree_model_get_path(mod, iter);
    s = gtk_tree_path_prev(path);
    if (s) {
	gtk_tree_model_get_iter(mod, iter, path);
    }
    gtk_tree_path_free(path);
    return s;
}

static void add_to_selection_count (GtkTreeModel *model, GtkTreePath *path,
				    GtkTreeIter *iter, int *count)
{
    *count += 1;
}

static void add_to_selection_list (GtkTreeModel *model, GtkTreePath *path,
				   GtkTreeIter *iter, int *list)
{
    gchar *varnum = NULL;

    gtk_tree_model_get(model, iter, 0, &varnum, -1);
    list[0] += 1;
    list[list[0]] = atoi(varnum);
    g_free(varnum);
}

int *main_window_selection_as_list (void) 
{
    GtkTreeSelection *select;
    int *list = NULL;
    int scount = 0;

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_selected_foreach(select, 
					(GtkTreeSelectionForeachFunc) 
					add_to_selection_count,
					&scount); 

    if (scount > 0) {
	list = gretl_list_new(scount);
    }

    if (list != NULL) {
	list[0] = 0;
	gtk_tree_selection_selected_foreach(select, 
					    (GtkTreeSelectionForeachFunc) 
					    add_to_selection_list,
					    list); 
    }

    return list;
}
