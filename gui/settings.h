#ifndef SETTINGS_H
#define SETTINGS_H

void set_rcfile (void);

void write_rc (void);

void mkfilelist (int filetype, const char *newfile);

void delete_from_filelist (int filetype, const char *fname);

void add_files_to_menu (int filetype);

void options_dialog (gpointer data);

void font_selector (void);

void load_fixed_font (void);

char *endbit (char *dest, char *src, int addscore);

void get_default_dir (char *s);

void filesel_set_path_callback (const char *setting, char *strvar);

#ifdef HAVE_TRAMO
void set_tramo_ok (int set);
#endif

#ifdef HAVE_X12A
void set_x12a_ok (int set);
#endif

#endif /* SETTINGS_H */
