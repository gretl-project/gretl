/* callbacks.h */

#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <gtk/gtk.h>

void listbox_select_row (GtkTreeSelection *selection, gpointer data);

gint listbox_double_click (GtkWidget *widget, GdkEventButton *event,
			   windata_t *win);

gboolean listbox_drag (GtkWidget *widget, GdkEventMotion *event,
		       gpointer data);
 
void open_data (gpointer data, guint dir_code, GtkWidget *widget);

void open_script (gpointer data, guint code, GtkWidget *widget);

void file_save (windata_t *vwin, guint file_code, GtkWidget *widget);

void dummy_call (void);

void print_report (gpointer data, guint unused, GtkWidget *widget);

void edit_header (gpointer data, guint save, GtkWidget *widget);

void fit_resid_callback (gpointer data, guint code, GtkWidget *widget);

void model_stat_callback (gpointer data, guint which, GtkWidget *widget);

void model_callback (gpointer data, guint model_code, GtkWidget *widget);

void selector_callback (gpointer data, guint action, GtkWidget *widget);

void gretl_callback (gpointer data, guint action, GtkWidget *widget);

void model_genr_callback (gpointer data, guint u, GtkWidget *widget);

void file_save_callback (GtkWidget *w, windata_t *vwin);

void newdata_callback (gpointer data, guint pd_code, GtkWidget *widget);

void xcorrgm_callback (gpointer p, guint v, GtkWidget *w);

void do_nistcheck (gpointer p, guint u, GtkWidget *w);

#if defined (ENABLE_MAILER) && !defined(G_OS_WIN32)
void send_file (char *fullname);
#endif

#endif /* CALLBACKS_H */
