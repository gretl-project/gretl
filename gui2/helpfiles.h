#ifndef HELPFILES_H
#define HELPFILES_H

void helpfile_init (void);

void context_help (GtkWidget *widget, gpointer data);

void do_gui_help (gpointer data, guint pos, GtkWidget *widget);

void do_script_help (gpointer data, guint pos, GtkWidget *widget);

gint edit_script_help (GtkWidget *widget, GdkEventButton *b,
			      windata_t *vwin);

void datafile_find (GtkWidget *widget, gpointer data);

void menu_find (gpointer data, guint dbfind, GtkWidget *widget);

void text_find_callback (GtkWidget *w, gpointer data);

int new_style_gui_help (FILE *fp);

#ifdef OLD_GTK
void colorize_tooltips (GtkTooltips *tip);
#endif

void gretl_tooltips_init (void);

void gretl_tooltips_add (GtkWidget *w, const gchar *str);

GtkItemFactoryEntry *get_help_menu_items (int code);

#endif /* HELPFILES_H */
