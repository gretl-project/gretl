#ifndef GUI_UTILS_H
#define GUI_UTILS_H

/* functions follow */

#if defined(G_OS_WIN32) || defined (USE_GNOME)
void window_print (windata_t *mydata, guint u, GtkWidget *widget);
#endif

#ifdef ENABLE_NLS
gchar *menu_translate (const gchar *path, gpointer p);
#endif

int getbufline (char *buf, char *line, int init);

void flip (GtkItemFactory *ifac, const char *path, gboolean s);

int copyfile (const char *src, const char *dest);

int isdir (const char *path);

void append_dir (char *fname, const char *dir);

void delete_widget (GtkWidget *widget, gpointer data);

void *mymalloc (size_t size); 

void *myrealloc (void *ptr, size_t size);

void mark_dataset_as_modified (void);

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
			char *title, int role,
			GtkItemFactoryEntry menu_items[]);

windata_t *view_file (char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role, 
		      GtkItemFactoryEntry menu_items[]);

windata_t *edit_buffer (char **pbuf, int hsize, int vsize, 
			char *title, int role);

int view_model (PRN *prn, MODEL *pmod, int hsize, int vsize, 
		char *title);

void errbox (const char *msg);

void infobox (const char *msg);

int validate_varname (const char *varname);

void text_copy (gpointer data, guint how, GtkWidget *widget);

void text_paste (windata_t *mydata, guint u, GtkWidget *widget);

void text_undo (windata_t *mydata, guint u, GtkWidget *widget);

gint popup_menu_handler (GtkWidget *widget, GdkEvent *event,
			 gpointer data);

void add_popup_item (const gchar *label, GtkWidget *menu,
		     GCallback callback, gpointer data);

void get_stats_table (void);

int gui_open_plugin (const char *plugin, void **handle);

void text_set_cursor (GtkWidget *w, GdkCursorType cspec);

gint get_char_width (GtkWidget *widget);

gchar *textview_get_text (GtkTextView *view);

int build_path (const char *dir, const char *fname, char *path, 
		const char *ext);

int prn_to_clipboard (PRN *prn, int copycode);

#define standard_button(f) gtk_button_new_from_stock(f)

int get_worksheet_data (const char *fname, int datatype, int append);

#endif /* GUI_UTILS_H */
