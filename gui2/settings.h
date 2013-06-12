#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef G_OS_WIN32
int read_win32_config (int debug);
void set_use_wimp (int s);
int get_use_wimp (void);
#else
int gretl_config_init (void);
#endif

#ifdef HAVE_TRAMO
int get_tramo_ok (void);
#endif

#ifdef HAVE_X12A
int get_x12a_ok (void);
#endif

#if defined(MAC_NATIVE) && defined(PKGBUILD)
void set_up_mac_look (void);
#endif

void set_gretl_startdir (void);

int using_hc_by_default (void);

int get_manpref (void);

int autoicon_on (void);

int use_tabbed_editor (void);

int use_tabbed_model_viewer (void);

int session_prompt_on (void);

void set_session_prompt (int val);

int get_keep_folder (void);

int write_rc (void);

void dump_rc (void);

void force_english_help (void);

int options_dialog (int page, const char *varname, GtkWidget *parent);

void font_selector (GtkAction *action);

void set_fixed_font (const char *fontname);

void update_persistent_graph_colors (void);

void update_persistent_graph_font (void);

void set_app_font (const char *fontname);

const char *get_app_fontname (void);

const char *get_fixed_fontname (void);

void get_default_dir (char *s, int action);

void working_dir_dialog (void);

int gui_set_working_dir (char *dirname);

void set_working_dir_callback (GtkWidget *w, char *path);

void set_path_callback (char *setvar, char *setting);

void set_datapage (const char *str);

void set_scriptpage (const char *str);

const char *get_datapage (void);

const char *get_scriptpage (void);

const char *get_default_hc_string (int ci);

int check_for_prog (const char *prog);

#endif /* SETTINGS_H */
