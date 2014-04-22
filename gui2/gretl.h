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

/* gretl.h, main header file for gretl gui */

#ifndef GRETL_H
#define GRETL_H

#ifdef _WIN32
# include "winconfig.h"
#else
# include "config.h"
#endif

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "dbread.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define GNULL (gconstpointer) NULL

/* fixups for GTK3 */

#if GTK_MAJOR_VERSION >= 3
#include <gdk/gdkkeysyms-compat.h>
# define gtk_combo_box_entry_new_text gtk_combo_box_text_new_with_entry
#endif

/* remedial macros for old GTK installations */

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# define gtk_widget_is_sensitive(w) GTK_WIDGET_IS_SENSITIVE(w)
# define gtk_widget_has_focus(w) GTK_WIDGET_HAS_FOCUS(w)
# define gtk_widget_set_can_default(w,s) GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT)
# define gtk_widget_set_can_focus(w,s) GTK_WIDGET_SET_FLAGS(w, GTK_CAN_FOCUS)
#endif

/* remedial macro for Mac */
#ifdef OS_OSX
# define right_click(e) (e->button == 3 || \
			 (e->button == 1 && \
			  (e->state & GDK_CONTROL_MASK)))
# define cmd_key(e) (e->state & GDK_META_MASK)
# define alt_w_key 0x1002211
# define alt_x_key 0x1002248
#else
# define right_click(e) (e->button == 3)
# define cmd_key(e) (0)
#endif

#include "gretltypes.h"
#include "gui_utils.h"
#include "callbacks.h"
#include "dialogs.h"
#include "library.h"
#include "settings.h"
#include "helpfiles.h"
#include "gui_recode.h"

#include "locale.h"

#define GRETL_BUFSIZE 8192
#define MAXSTR FILENAME_MAX
#define CLIPTEMP "cliptmp.txt"

#define SCRIPT_WIDTH 78
#define SCRIPT_HEIGHT 370
#define MODEL_WIDTH 72
#define MODEL_HEIGHT 420 

/* basic global program vars */
extern DATASET *dataset;
extern MODEL *model;

/* global counters */
extern int orig_vars;

/* global state variables */
extern int data_status;
extern float gui_scale;

/* global filenames */
extern char datafile[MAXLEN];
extern char scriptfile[MAXLEN];
extern char tryfile[MAXLEN];

/* global option-related vars */
extern int winsize;
extern int main_x;
extern int main_y;
extern int mainwin_width;
extern int mainwin_height;
extern int ox_support;

#if !defined(G_OS_WIN32) && !defined(OS_OSX)
extern char viewps[MAXSTR];
extern char viewpdf[MAXSTR];
extern char Browser[MAXSTR];
#endif

extern char calculator[MAXSTR];
extern char latex[MAXSTR];
extern char viewdvi[MAXSTR];
extern char Rcommand[MAXSTR];

#if defined(HAVE_AUDIO) && !defined(G_OS_WIN32)
extern char midiplayer[MAXSTR];
#endif

/* global GUI equipment */
extern windata_t *mdata;
extern GtkTargetEntry gretl_drag_targets[];
extern PangoFontDescription *fixed_font;

#include "gretl_enums.h"

/* functions follow */

#ifndef WIN32
int gretl_fork (const char *progvar, const char *arg);
void set_wm_icon (GtkWidget *w);
#endif

int mdata_selection_count (void);
int mdata_active_var (void);
void populate_varlist (void);
void clear_varlist (GtkWidget *widget);
void mdata_select_last_var (void);
int gui_restore_sample (DATASET *dset);
void make_list_from_main (void);
void do_stop_script (GtkWidget *w, windata_t *vwin);
int is_control_key (guint k);
gboolean real_open_tryfile (void);
gchar *user_friendly_menu_path (const char *mpath,
				gboolean modelwin);

#ifdef MAC_INTEGRATION
gint mac_hide_unhide (GdkEventKey *event);
#endif

/* functions defined in files other than gretl.c */
void about_dialog (void);

#endif /* GRETL_H */
