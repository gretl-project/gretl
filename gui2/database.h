#ifndef GRETL_DATABASE_H
#define GRETL_DATABASE_H

void open_db_list (GtkWidget *w, gpointer data);
void open_named_db_list (char *dbname);
void open_remote_db_list (GtkWidget *w, gpointer data);
void open_named_remote_db_list (char *dbname);
void grab_remote_db (GtkWidget *w, gpointer data);
gint populate_dbfilelist (windata_t *ddata);
void display_db_error (windata_t *dbwin, char *buf);
void do_compact_data_set (void);
void do_expand_data_set (void);
gchar *get_db_description (const char *binname);
int write_db_description (const char *binname, const char *descrip);

#endif
