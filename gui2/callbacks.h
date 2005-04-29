/* callbacks.h */

#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <gtk/gtk.h>

/* functions follow */

#if GTK_MAJOR_VERSION == 1
void selectrow (GtkCList *clist, gint row, gint column, 
                GdkEventButton *event, gpointer data);
void unselectrow (GtkCList *clist, gint row, gint column, 
		  GdkEventButton *event, gpointer data);
#else

void listbox_select_row (GtkTreeSelection *selection, gpointer data);

gint listbox_double_click (GtkWidget *widget, GdkEventButton *event,
			   windata_t *win);

gboolean listbox_drag (GtkWidget *widget, GdkEventMotion *event,
		       gpointer data);
#endif
 
void open_data (gpointer data, guint dir_code, GtkWidget *widget);

void open_script (gpointer data, guint code, GtkWidget *widget);

void file_save (gpointer data, guint file_code, GtkWidget *widget);

void dummy_call (void);

void print_report (gpointer data, guint unused, GtkWidget *widget);

void edit_header (gpointer data, guint save, GtkWidget *widget);

void fit_resid_callback (gpointer data, guint code, GtkWidget *widget);

void var_resid_callback (gpointer data, guint eqnum, GtkWidget *widget);

void model_stat_callback (gpointer data, guint which, GtkWidget *widget);

void model_callback (gpointer data, guint model_code, GtkWidget *widget);

#ifdef ENABLE_GMP
void mp_ols_callback (gpointer data, guint model_code, GtkWidget *widget);
#endif

void selector_callback (gpointer data, guint action, GtkWidget *widget);

void gretl_callback (gpointer data, guint action, GtkWidget *widget);

void model_genr_callback (gpointer data, guint u, GtkWidget *widget);

void text_copy_callback (GtkWidget *w, gpointer data);

void text_paste_callback (GtkWidget *w, gpointer data);

void text_replace_callback (GtkWidget *w, gpointer data);

void text_undo_callback (GtkWidget *w, gpointer data);

void run_script_callback (GtkWidget *w, gpointer data);

void gp_send_callback (GtkWidget *w, gpointer data);

void file_save_callback (GtkWidget *w, gpointer data);

void add_random_callback (gpointer data, guint code, GtkWidget *widget);

void newdata_callback (gpointer data, guint pd_code, GtkWidget *widget);

#if 0
void start_panel_callback (gpointer data, guint u, GtkWidget *widget);
#endif

void do_nistcheck (gpointer p, guint u, GtkWidget *w);

#endif /* CALLBACKS_H */
