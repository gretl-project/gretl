#ifndef GRETL_DATABASE_H
#define GRETL_DATABASE_H

void open_db_list (GtkWidget *w, gpointer data);
void open_remote_db_list (GtkWidget *w, gpointer data);
void grab_remote_db (GtkWidget *w, gpointer data);
gint populate_dbfilelist (windata_t *ddata);

#endif
