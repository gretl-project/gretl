/* callbacks.h */

#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <gtk/gtk.h>

/* functions follow */
 
void selectrow (GtkCList *clist, gint row, gint column, 
	        GdkEventButton *event, gpointer data) ;

void open_data (gpointer data, guint dir_code, GtkWidget *widget);

void open_script (gpointer data, guint code, GtkWidget *widget);

void file_save (gpointer data, guint file_code, GtkWidget *widget);

void dummy_call (void);

void edit_header (gpointer data, guint save, GtkWidget *widgetvoid);

void fit_resid_callback (gpointer data, guint code, GtkWidget *widget);

void model_stat_callback (gpointer data, guint which, GtkWidget *widget);

void model_callback (gpointer data, guint model_code, GtkWidget *widget);

void gretl_callback (gpointer data, guint action, GtkWidget *widget);

void model_test_callback (gpointer data, guint action, GtkWidget *widget);

void text_copy_callback (GtkWidget *w, gpointer data);

void text_paste_callback (GtkWidget *w, gpointer data);

void text_replace_callback (GtkWidget *w, gpointer data);

void text_undo_callback (GtkWidget *w, gpointer data);

void run_script_callback (GtkWidget *w, gpointer data);

void file_save_callback (GtkWidget *w, gpointer data);

#endif /* CALLBACKS_H */
