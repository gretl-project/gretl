/* dialogs.h */

#ifndef DIALOGS_H
#define DIALOGS_H

enum {
    GRETL_YES,
    GRETL_NO,
    GRETL_CANCEL
} buttons;

/* functions follow */

int make_default_storelist (void);

void random_dialog (gpointer data, guint uni, GtkWidget *widget);

void newdata_dialog (gpointer data, guint pd_code, GtkWidget *widget);

void start_panel_dialog (gpointer data, guint pd_code, GtkWidget *widget);

void addvars_dialog (gpointer data, guint add_code, GtkWidget *widget);

void edit_dialog (char *diagtxt, char *infotxt, char *deftext, 
		  char *oktxt, void (*okfunc)(), void *okptr,
		  char *canceltxt, guint hlpcode, guint varclick);

void about_dialog (gpointer data);

gint yes_no_dialog (char *title, char *msg, int cancel);

void destroy_dialog_data (GtkWidget *w, gpointer data);

gint exit_check (GtkWidget *widget, GdkEvent *event, gpointer data);

void menu_exit_check (GtkWidget *w, gpointer data);

void delimiter_dialog (void);

void copy_format_dialog (windata_t *vwin);

void varinfo_dialog (int varnum);

#endif /* DIALOGS_H */
