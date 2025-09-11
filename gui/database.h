#ifndef GRETL_DATABASE_H
#define GRETL_DATABASE_H

void open_db_index (GtkWidget *w, gpointer data);

gboolean open_named_db_index (char *dbname);

void open_remote_db_index (GtkWidget *w, gpointer data);

gboolean open_named_remote_db_index (char *dbname);

void open_dbnomics_provider (GtkWidget *w, gpointer data);

void open_dbnomics_category (GtkWidget *w, gpointer data);

void open_dbnomics_dataset (GtkWidget *w, gpointer data);

void show_dbnomics_dimensions (GtkWidget *w, gpointer data);

void open_dbnomics_series (GtkWidget *w, gpointer data);

void dbnomics_specific_series (GtkAction *action, gpointer data);

void dbnomics_pager_call (GtkWidget *w, windata_t *vwin);

void install_file_from_server (GtkWidget *w, windata_t *vwin);

void drag_file_from_server (guint info);

void pkg_info_from_server (GtkWidget *w, windata_t *vwin);

void maybe_update_pkgview (const char *filename,
			   const char *pkgname,
			   int zipfile,
			   GtkWidget *parent);

gint populate_dbfilelist (windata_t *vwin, int *pndb);

void set_db_dir_callback (windata_t *vwin, char *path);

gint populate_remote_db_list (windata_t *vwin);

gint populate_dbnomics_provider_list (windata_t *vwin);

gint populate_dbnomics_category_list (windata_t *vwin,
                                      gchar *path,
                                      void *data);

gint populate_dbnomics_dataset_list (windata_t *vwin,
                                     gchar *path,
                                     void *data);

gint populate_dbnomics_series_list (windata_t *vwin,
                                    gchar *path);

gint populate_remote_func_list (windata_t *win, int filter);

gint populate_addons_list (windata_t *vwin);

gint populate_remote_data_pkg_list (windata_t *vwin);

void display_db_series (windata_t *vwin);

void drag_import_db_series (void);

void do_compact_dataset (void);

void do_expand_dataset (void);

gchar *get_db_description (const char *binname);

int write_db_description (const char *binname, const char *descrip);

void show_network_error (windata_t *vwin);

void open_rats_window (char *fname);

void open_bn7_window (char *fname);

void sync_db_windows (void);

int add_dbnomics_data (windata_t *vwin);

int show_dbnomics_data (windata_t *vwin, int plot);

void dbnomics_search (gchar *key, windata_t *vwin, int stype);

void maybe_fill_dbn_finder (GtkWidget *entry);

#endif
