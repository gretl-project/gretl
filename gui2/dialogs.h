/* dialogs.h */

#ifndef DIALOGS_H
#define DIALOGS_H

enum {
    GRETL_YES,
    GRETL_NO,
    GRETL_CANCEL,
    HELP_BUTTON
} buttons;

/* functions follow */

int make_default_storelist (void);

int storevars_dialog (int export); 

void random_dialog (gpointer data, guint uni, GtkWidget *widget);

void newdata_dialog (gpointer data, guint pd_code, GtkWidget *widget);

void start_panel_dialog (gpointer data, guint pd_code, GtkWidget *widget);

void edit_dialog (const char *diagtxt, const char *infotxt, const char *deftext, 
		  void (*okfunc)(), void *okptr,
		  guint hlpcode, guint varclick);

void about_dialog (gpointer data);

gint yes_no_dialog (char *title, char *msg, int cancel);

void destroy_dialog_data (GtkWidget *w, gpointer data);

gint exit_check (GtkWidget *widget, GdkEvent *event, gpointer data);

void menu_exit_check (GtkWidget *w, gpointer data);

void delimiter_dialog (void);

void copy_format_dialog (windata_t *vwin);

void varinfo_dialog (int varnum, int full);

#endif /* DIALOGS_H */
