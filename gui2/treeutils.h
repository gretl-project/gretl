#ifndef TREEUTILS_H
#define TREEUTILS_H

GtkWidget *list_box_create (windata_t *win, GtkBox *box, 
			    gint ncols, const char *titles[]);

void tree_view_get_string (GtkTreeView *view, int row, int col, gchar **val);

int tree_path_get_row_number (GtkTreePath *path); 

gboolean main_varclick (GtkWidget *widget, GdkEventButton *event,
			windata_t *win);

#endif /* TREEUTILS_H */
