#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef G_OS_WIN32
void read_rc (void);
#endif

void set_rcfile (void);

void write_rc (void);

void mkfilelist (int filetype, const char *newfile);

void delete_from_filelist (int filetype, const char *fname);

void add_files_to_menu (int filetype);

void options_dialog (gpointer data);

void font_selector (gpointer data, guint fixed, GtkWidget *widget);

void set_fixed_font (void);

void set_app_font (const char *fontname);

char *endbit (char *dest, char *src, int addscore);

void get_default_dir (char *s);

#endif /* SETTINGS_H */
