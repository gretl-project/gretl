#ifndef GRETL_DATABASE_H
#define GRETL_DATABASE_H

void open_db_index (GtkWidget *w, gpointer data);

void open_named_db_index (char *dbname);

void open_remote_db_index (GtkWidget *w, gpointer data);

void open_named_remote_db_index (char *dbname);

void install_file_from_server (GtkWidget *w, gpointer data);

void file_info_from_server (GtkWidget *w, gpointer data);

gint populate_dbfilelist (windata_t *ddata);

gint populate_remote_db_list (windata_t *vwin);

gint populate_remote_func_list (windata_t *win);

void gui_get_db_series (gpointer p, guint action, GtkWidget *w);

void import_db_series (windata_t *vwin);

void do_compact_data_set (void);

void do_expand_data_set (void);

gchar *get_db_description (const char *binname);

int write_db_description (const char *binname, const char *descrip);

void show_network_error (windata_t *vwin);

void open_rats_window (char *fname);

void open_bn7_window (char *fname);

#endif
