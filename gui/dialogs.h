/* dialogs.h */

#ifndef DIALOGS_H
#define DIALOGS_H

/* functions follow */

int make_default_storelist (void);

int storevars_dialog (int export); 

void random_dialog (gpointer data, guint uni, GtkWidget *widget);

void newdata_dialog (gpointer data, guint pd_code, GtkWidget *widget);

void start_panel_dialog (gpointer data, guint pd_code, GtkWidget *widget);

void addvars_dialog (gpointer data, guint add_code, GtkWidget *widget);

void graph_dialog (gpointer data, guint controls, GtkWidget *widget); 

void edit_dialog (char *diagtxt, char *infotxt, char *deftext, 
		  int edit_shown,
		  char *oktxt, void (*okfunc)(), void *okptr,
		  char *canceltxt, void (*cancelfunc)(), 
		  void *cancelptr, guint hlpcode, guint varclick);

void about_dialog (gpointer data);

gint yes_no_dialog (char *title, char *msg, int cancel);

void destroy_dialog_data (GtkWidget *w, gpointer data);

gint exit_check (GtkWidget *widget, GdkEvent *event, gpointer data);

void menu_exit_check (GtkWidget *w, gpointer data);

#endif /* DIALOGS_H */
