#ifndef TREEUTILS_H
#define TREEUTILS_H

int tree_selection_count (GtkTreeSelection *select, int *vnum);

int vwin_selection_count (windata_t *vwin, int *row);

void vwin_add_list_box (windata_t *vwin, GtkBox *box, 
			int ncols, gboolean hidden_col,
			GType *types, const char **titles,
			int tree);

void presort_treelist (windata_t *vwin);

void tree_view_get_bool (GtkTreeView *view, int row, int col, gboolean *val);

void tree_view_get_string (GtkTreeView *view, int row, int col, gchar **val);

void list_store_set_string (GtkTreeView *view, int row, int col, const gchar *val);

void tree_store_set_string (GtkTreeView *view, int row, int col, const gchar *val);

void tree_view_get_int (GtkTreeView *view, int row, int col, int *val);

int tree_path_get_row_number (GtkTreePath *path); 

void tree_model_get_iter_last (GtkTreeModel *mod, GtkTreeIter *iter);

gboolean tree_model_iter_prev (GtkTreeModel *mod, GtkTreeIter *iter);

gboolean main_varclick (GtkWidget *widget, GdkEventButton *event,
			windata_t *win);

int *main_window_selection_as_list (void); 

void set_main_colheads_clickable (gboolean s);

#endif /* TREEUTILS_H */
