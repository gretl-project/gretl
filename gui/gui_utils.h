#ifndef GUI_UTILS_H
#define GUI_UTILS_H

/* functions follow */

#ifdef USE_GNOME
void window_print (windata_t *mydata, guint u, GtkWidget *widget);
#endif

#ifdef ENABLE_NLS
gchar *menu_translate (const gchar *path, gpointer p);
#endif

int getbufline (char *buf, char *line, int init);

void flip (GtkItemFactory *ifac, char *path, gboolean s);

int copyfile (const char *src, const char *dest);

int prn_to_clipboard (PRN *prn, int copycode);

int isdir (const char *path);

void delete_model (GtkWidget *widget, gpointer data);

void delete_widget (GtkWidget *widget, gpointer data);

void catch_view_key (GtkWidget *w, GdkEventKey *key);

void *mymalloc (size_t size); 

void *myrealloc (void *ptr, size_t size);

void mark_dataset_as_modified (void);

void clear_data (void);

void register_data (const char *fname, const char *user_fname,
		    int record);

void do_open_data (GtkWidget *w, gpointer data, int code);

void verify_open_data (gpointer userdata, int code);

void verify_open_session (gpointer userdata);

void save_session (char *fname);

void close_window (gpointer data, guint win_code, GtkWidget *widget);

void windata_init (windata_t *mydata);

void free_windata (GtkWidget *w, gpointer data);

windata_t *view_buffer (PRN *prn, int hsize, int vsize, 
			const char *title, int role,
			gpointer data);

windata_t *view_file (char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role); 

windata_t *edit_buffer (char **pbuf, int hsize, int vsize, 
			char *title, int role);

int view_model (PRN *prn, MODEL *pmod, int hsize, int vsize, 
		char *title);

void file_view_set_editable (windata_t *vwin);

void setup_column (GtkWidget *listbox, int column, int width);

void errbox (const char *msg);

void infobox (const char *msg);

int validate_varname (const char *varname);

void text_copy (gpointer data, guint how, GtkWidget *widget);

void text_paste (windata_t *mydata, guint u, GtkWidget *widget);

void text_undo (windata_t *mydata, guint u, GtkWidget *widget);

gint popup_menu_handler (GtkWidget *widget, GdkEvent *event, gpointer data);

void add_popup_item (gchar *label, GtkWidget *menu,
		     GtkSignalFunc func, gpointer data);

void get_stats_table (void);

void *gui_get_plugin_function (const char *funcname, 
			       void **phandle);

int build_path (const char *dir, const char *fname, char *path, const char *ext);

int get_worksheet_data (const char *fname, int datatype, int append);

char *double_underscores (char *targ, const char *src);

#endif /* GUI_UTILS_H */
