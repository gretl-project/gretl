/* dialogs.h */

#ifndef DIALOGS_H
#define DIALOGS_H

enum {
    GRETL_YES,
    GRETL_NO,
    GRETL_CANCEL
} buttons;

/* functions follow */

void errbox (const char *msg);

void infobox (const char *msg);

int make_default_storelist (void);

void addvars_dialog (gpointer data, guint add_code, GtkWidget *widget);

void edit_dialog (char *diagtxt, char *infotxt, char *deftext, 
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

void set_sample_dialog (gpointer p, guint u, GtkWidget *w);

void arma_options_dialog (gpointer p, guint u, GtkWidget *w);

GtkWidget *standard_button (int code);

#endif /* DIALOGS_H */
