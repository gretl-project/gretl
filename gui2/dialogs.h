/* dialogs.h */

#ifndef DIALOGS_H
#define DIALOGS_H

typedef struct dialog_t_ dialog_t;

enum {
    GRETL_YES,
    GRETL_NO,
    GRETL_CANCEL,
    HELP_BUTTON
} buttons;

extern GtkWidget *open_dialog;

/* functions follow */

void errbox (const char *msg);

void infobox (const char *msg);

gint yes_no_dialog (char *title, char *msg, int cancel);

int make_default_storelist (void);

void edit_dialog (const char *diagtxt, const char *infotxt, const char *deftext, 
		  void (*okfunc)(), void *okptr,
		  guint hlpcode, guint varclick);

const gchar *dialog_data_get_text (dialog_t *ddata);

gchar *dialog_data_special_get_text (dialog_t *ddata);

int dialog_data_get_action (const dialog_t *ddata);

gretlopt dialog_data_get_opt (const dialog_t *ddata);

gpointer dialog_data_get_data (dialog_t *ddata);

GtkWidget *dialog_data_get_vbox (dialog_t *ddata);

void close_dialog (dialog_t *ddata);

gint exit_check (GtkWidget *widget, GdkEvent *event, gpointer data);

void menu_exit_check (GtkWidget *w, gpointer data);

void delimiter_dialog (void);

void copy_format_dialog (windata_t *vwin, int multicopy);

void varinfo_dialog (int varnum, int full);

void sample_range_dialog (gpointer p, guint u, GtkWidget *w);

void arma_options_dialog (gpointer p, guint u, GtkWidget *w);

void panel_structure_dialog (DATAINFO *pdinfo, GtkWidget *w);

void data_compact_dialog (GtkWidget *w, int spd, int *target_pd, 
			  int *mon_start, gint *compact_method);


#ifdef OLD_GTK
GtkWidget *standard_button (int code);
#endif

#endif /* DIALOGS_H */
