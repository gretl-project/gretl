/* callbacks.h */

#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <gtk/gtk.h>

void listbox_select_row (GtkTreeSelection *selection, gpointer data);

gint listbox_double_click (GtkWidget *widget, GdkEventButton *event,
			   windata_t *win);

gboolean listbox_drag (GtkWidget *widget, GdkEventMotion *event,
		       gpointer data);
 
void open_data (GtkAction *action);

void file_save (windata_t *vwin, int ci);

void fsave_callback (GtkAction *action, gpointer p);

void dummy_call (void);

void fit_resid_callback (GtkAction *action, gpointer data);

void model_stat_callback (GtkAction *action, gpointer data);

void model_callback (GtkAction *action, gpointer data);

void selector_callback (GtkAction *action, gpointer data);

void gretl_callback (GtkAction *action, gpointer data);

void model_genr_callback (GtkAction *action, gpointer data);

void revise_nl_model (MODEL *pmod, GtkWidget *parent);

void revise_system_model (void *ptr, GtkWidget *parent);

void newdata_callback (void);

void xcorrgm_callback (void);

void do_nistcheck (GtkAction *action);

void genr_callback (void);

void minibuf_callback (void);

void menu_boxplot_callback (int varnum);

void boxplot_callback (void);

#if defined (ENABLE_MAILER) && !defined(G_OS_WIN32)
void send_file (char *fullname);
#endif

#endif /* CALLBACKS_H */
