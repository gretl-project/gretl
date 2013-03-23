#ifndef GRETL_DATABASE_H
#define GRETL_DATABASE_H

void open_db_index (GtkWidget *w, gpointer data);

gboolean open_named_db_index (char *dbname);

void open_remote_db_index (GtkWidget *w, gpointer data);

gboolean open_named_remote_db_index (char *dbname);

void install_file_from_server (GtkWidget *w, windata_t *vwin);

void pkg_info_from_server (GtkWidget *w, windata_t *vwin);

int unzip_package_file (const char *zipname, const char *path);

gint populate_dbfilelist (windata_t *vwin, int *pndb);

void set_db_dir_callback (windata_t *vwin, char *path);

gint populate_remote_db_list (windata_t *vwin);

gint populate_remote_func_list (windata_t *win);

gint populate_remote_addons_list (windata_t *vwin);

gint populate_remote_data_pkg_list (windata_t *vwin);

void display_db_series (windata_t *vwin);

void import_db_series (windata_t *vwin);

void do_compact_data_set (void);

void do_expand_data_set (void);

gchar *get_db_description (const char *binname);

int write_db_description (const char *binname, const char *descrip);

void show_network_error (windata_t *vwin);

void open_rats_window (char *fname);

void open_bn7_window (char *fname);

void sync_db_windows (void);

#endif
