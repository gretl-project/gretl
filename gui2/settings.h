#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef G_OS_WIN32
void read_rc (void);
#else
void gretl_config_init (void);
#endif

void set_gretl_startdir (void);

int using_hc_by_default (void);

int get_manpref (void);

void write_rc (void);

void dump_rc (void);

void force_english_help (void);

int options_dialog (int page, const char *varname);

void font_selector (gpointer data, guint which, GtkWidget *widget);

void set_fixed_font (void);

#ifndef USE_GNOME

void set_app_font (const char *fontname);

const char *get_app_fontname (void);

#endif

void graph_color_selector (GtkWidget *w, gpointer p);

GtkWidget *color_patch_button (int cnum);

void color_patch_button_reset (GtkWidget *button, int cnum);

void get_default_dir (char *s, int action);

int gui_set_working_dir (char *dirname);

void set_working_dir_callback (GtkWidget *w, char *path);

void working_dir_dialog (void);

void set_path_callback (char *setvar, char *setting);

void set_datapage (const char *str);

void set_scriptpage (const char *str);

const char *get_datapage (void);

const char *get_scriptpage (void);

int check_for_prog (const char *prog);

#endif /* SETTINGS_H */
