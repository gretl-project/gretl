/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef G_OS_WIN32
int read_win32_config (int debug, int ignore_rc);
#else
int gretl_config_init (int ignore_rc);
#endif

#ifdef HAVE_TRAMO
int get_tramo_ok (void);
#endif

#ifdef HAVE_X12A
int get_x12a_ok (void);
#endif

#if defined(__APPLE__) && defined(HAVE_MAC_THEMES)
void set_up_mac_look (void);
#endif

#if defined(G_OS_WIN32)
void set_up_windows_look (void);
# if GTK_MAJOR_VERSION < 3
void set_wimp_preferred (int s);
# endif
#endif

void set_gretl_startdir (void);

int using_hc_by_default (void);

int get_manpref (void);

int autoicon_on (void);

int get_icon_sizing (void);

int dark_theme_active (void);

const char *blue_for_text (void);

int use_tabbed_editor (void);

int use_tabbed_model_viewer (void);

int session_prompt_on (void);

void set_session_prompt (int val);

int get_keep_folder (void);

void set_script_output_policy (int p, windata_t *vwin);

int get_script_output_policy (void);

int write_rc (gretlopt opt);

void sync_path_from_lib (const char *path_id);

void dump_rc (void);

void force_english_help (void);

int preferences_dialog (int page, const char *varname, GtkWidget *parent);

int console_prefs_dialog (GtkWidget *caller);

void font_selector (GtkAction *action);

void font_scale_selector (GtkAction *action);

void set_fixed_font (const char *fontname, int remember);

void set_app_font (const char *fontname, int remember);

const char *get_app_fontname (void);

const char *get_fixed_fontname (void);

void update_persistent_graph_colors (void);

void update_persistent_graph_font (void);

void get_default_dir_for_action (char *s, int action);

#ifndef GRETL_EDIT

void workdir_dialog0 (void);

void workdir_dialog1 (void);

int gui_set_working_dir (char *dirname);

void set_working_dir_callback (GtkWidget *w, char *path);

double next_graph_scale (double s, int mod);

double min_graph_scale (void);

double max_graph_scale (void);

#endif

void set_path_callback (char *setvar, char *setting);

void set_datapage (const char *str);

void set_scriptpage (const char *str);

const char *get_datapage (void);

const char *get_scriptpage (void);

const char *get_default_hc_string (int ci);

int check_for_prog (const char *prog);

void get_model_table_prefs (int *colheads,
			    int *use_tstats,
			    int *do_pvals,
			    int *do_asts,
			    int *figs,
			    char *fmt);

void set_model_table_prefs (int colheads,
			    int use_tstats,
			    int do_pvals,
			    int do_asts,
			    int figs,
			    char fmt);

void set_author_mail (const char *s);

const char *get_author_mail (void);

const char *get_sourceview_style (void);

#endif /* SETTINGS_H */
