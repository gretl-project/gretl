#ifndef FILELISTS_H
#define FILELISTS_H

#define MAXRECENT 4

void mkfilelist (int filetype, char *newfile);

void init_fileptrs (void);

void initialize_file_lists (void);

void write_filename_to_list (int filetype, int i, char *fname);

void delete_from_filelist (int filetype, const char *fname);

void add_files_to_menus (void);

void save_file_lists (FILE *fp);

char *endbit (char *dest, char *src, int addscore);


#endif /* FILELISTS_H */
