#ifndef FILELISTS_H
#define FILELISTS_H

#define MAXRECENT 4

void mkfilelist (int filetype, char *newfile);

void init_fileptrs (void);

void initialize_file_lists (void);

void write_filename_to_list (int filetype, int i, char *fname);

void delete_from_filelist (int filetype, const char *fname);

void add_files_to_menus (void);

char *endbit (char *dest, char *src, int addscore);

#if defined(USE_GNOME) && !defined(OLD_GTK)
# include <gconf/gconf-client.h>
void save_file_lists (GConfClient *client);
#elif defined(USE_GNOME) || defined(G_OS_WIN32)
void save_file_lists (void);
#else
void save_file_lists (FILE *fp);
#endif


#endif /* FILELISTS_H */
