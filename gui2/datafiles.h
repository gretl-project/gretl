#ifndef GRETL_DATAFILES_H
#define GRETL_DATAFILES_H

void browser_open_data (GtkWidget *w, gpointer data);
void browser_open_ps (GtkWidget *w, gpointer data);
void destroy_file_collections (void);
gint populate_filelist (windata_t *fdata, gpointer p);
char *strip_extension (char *s);
void panel_structure_dialog (DATAINFO *pdinfo, GtkWidget *w);

#endif
