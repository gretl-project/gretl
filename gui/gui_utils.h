#ifndef GUI_UTILS_H
#define GUI_UTILS_H

/* functions follow */

#ifdef G_OS_WIN32
void read_rc (void);
#endif

void write_rc (void);

int getbufline (char *buf, char *line, int init);

void flip (GtkItemFactory *ifac, char *path, gboolean s);

void mkfilelist (int filetype, const char *newfile);

void add_files_to_menu (int filetype);

int copyfile (const char *src, const char *dest);

int isdir (const char *path);

void append_dir (char *fname, const char *dir);
 
void set_rcfile (void);

void delete_model (GtkWidget *widget, gpointer data);

void delete_widget (GtkWidget *widget, gpointer data);

void catch_key (GtkWidget *w, GdkEventKey *key);

void *mymalloc (size_t size); 

void *myrealloc (void *ptr, size_t size);

void clear_data (void);

void register_data (const char *fname);

void verify_open_data (gpointer userdata);

void verify_open_session (gpointer userdata);

void save_session (char *fname);

void helpwin (gpointer data, guint script, GtkWidget *widget);

void menu_find (gpointer data, guint dbfind, GtkWidget *widget);

void close_window (gpointer data, guint win_code, GtkWidget *widget);

void context_help (GtkWidget *widget, gpointer data);

void do_help (gpointer data, guint help_code, GtkWidget *widget);

void windata_init (windata_t *mydata);

void free_windata (GtkWidget *w, gpointer data);

int view_buffer (print_t *prn, int hsize, int vsize, char *title, 
		 GtkItemFactoryEntry menu_items[], int menusize);

int view_file (char *filename, int editable, int del_file, 
	       int hsize, int vsize, char *title, 
	       GtkItemFactoryEntry menu_items[], int menusize);

int view_model (void *ptr, MODEL *pmod, int hsize, int vsize, 
		char *title);

void setup_column (GtkWidget *listbox, int column, int width);

void errbox (const char *msg);

void infobox (const char *msg);

int validate_varname (const char *varname);

void options_dialog (gpointer data);

void font_selector (void);

void text_copy (gpointer data, guint all, GtkWidget *widget);

void make_menu_item (gchar *label, GtkWidget *menu,
		     GtkSignalFunc func, gpointer data);

void get_stats_table (void);

int open_plugin (const char *plugin, void **handle);

void *get_plugin_function (const char *funcname, void *handle);

void close_plugin (void *handle);

#endif /* GUI_UTILS_H */
