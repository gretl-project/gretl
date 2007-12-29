#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef G_OS_WIN32
void read_rc (void);
#else
void gretl_config_init (void);
#endif

void set_program_startdir (void);

int using_hc_by_default (void);

int get_manpref (void);

void write_rc (void);

void dump_rc (void);

void force_english_help (void);

int options_dialog (int page);

void options_dialog_callback (gpointer p, guint u, GtkWidget *w);

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

void gui_set_working_dir (char *dirname);

void set_path_callback (char *setvar, char *setting);

void set_datapage (const char *str);

void set_scriptpage (const char *str);

const char *get_datapage (void);

const char *get_scriptpage (void);

#ifndef G_OS_WIN32
void first_time_set_user_dir (void);
#endif

int check_for_prog (const char *prog);

#endif /* SETTINGS_H */
