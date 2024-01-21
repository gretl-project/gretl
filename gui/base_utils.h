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

#ifndef BASE_UTILS_H
#define BASE_UTILS_H

enum { ADD_COMMA, DROP_COMMA };

gchar *gretl_window_title (const char *s);

int copyfile (const char *src, const char *dest);

int gretl_file_get_contents (const gchar *fname,
			     gchar **contents,
			     gsize *size);

FILE *gretl_tempfile_open (char *fname);

int bufopen (PRN **pprn);

PRN *gui_prn_new (void);

void set_wait_cursor (GdkWindow **pcwin);

void unset_wait_cursor (GdkWindow *cwin);

void gretl_set_window_modal (GtkWidget *w);

void gretl_set_window_quasi_modal (GtkWidget *w);

void dummy_call (void);

void nomem (void);

void *mymalloc (size_t size);

void *myrealloc (void *ptr, size_t size);

const char *path_last_slash_const (const char *path);

char *gretl_basename (char *dest, const char *src,
		      int addscore);

char *double_underscores (char *targ, const char *src);

gchar *double_underscores_new (const char *src);

char *adjust_fontspec_string (char *targ, const char *src,
			      int mod);

int get_string_width (const gchar *str);

int font_has_symbol (PangoFontDescription *desc, int symbol);

const char *print_today (void);

void *gui_get_plugin_function (const char *funcname);

gboolean do_open_script (int action);

void do_new_script (int code, const char *buf,
		    const char *scriptname);

void do_run_script (GtkWidget *w, windata_t *vwin);

void run_script_silent (GtkWidget *w, windata_t *vwin);

void new_script_callback (GtkAction *action);

void start_R (const char *buf, int send_data, int interactive);

#ifndef G_OS_WIN32

int browser_open (const char *url);

int gretl_fork (const char *progvar, const char *arg,
		const char *opt);

#endif

#endif /* BASE_UTILS_H */
