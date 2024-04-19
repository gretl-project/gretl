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

#include "config.h"

#define FULL_XML_HEADERS
#include "libgretl.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#include <string.h>

#define GNULL (gconstpointer) NULL

/* fixups for GTK3 */

#if GTK_MAJOR_VERSION >= 3
# include <gdk/gdkkeysyms-compat.h>
# include <gtksourceview/gtksource.h>
# define gtk_combo_box_entry_new_text gtk_combo_box_text_new_with_entry
# define GTKSOURCEVIEW_VERSION GTK_SOURCE_MAJOR_VERSION
#else
# include <gtksourceview/gtksourceview.h>
# define GTKSOURCEVIEW_VERSION 2
#endif

/* convenience macros */
#define widget_get_int(w,s) GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), s))
#define widget_set_int(w,s,i) g_object_set_data(G_OBJECT(w), s, GINT_TO_POINTER(i))

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
#include "base_utils.h"
#include "viewers.h"
#include "callbacks.h"
#include "dialogs.h"
#include "library.h"
#include "settings.h"
#include "helpfiles.h"
#include "focus.h"
#include "tabwin.h"

#ifdef ENABLE_NLS
# include "locale.h"
#endif

#define GRETL_BUFSIZE 8192
#define MAXSTR FILENAME_MAX
#define CLIPTEMP_TXT "cliptmp.txt"
#define CLIPTEMP_GDT "cliptmp.gdt"

#define SCRIPT_WIDTH 78
#define SCRIPT_HEIGHT 420
#define MODEL_WIDTH 72
#define MODEL_HEIGHT 420

#ifndef GRETL_EDIT
/* basic global program vars */
extern DATASET *dataset;
extern MODEL *model;
#endif

/* global counters */
extern int orig_vars;

/* global state variables */
extern int data_status;
extern float gui_scale;

/* global filenames */
extern char datafile[MAXLEN];
extern char scriptfile[MAXLEN];

/* global option-related vars */
extern int winsize;
extern int main_x;
extern int main_y;
extern int mainwin_width;
extern int mainwin_height;
extern int swallow;

#if !defined(G_OS_WIN32) && !defined(OS_OSX)
extern char viewps[MAXSTR];
extern char viewpdf[MAXSTR];
extern char Browser[MAXSTR];
#endif

extern char calculator[MAXSTR];
extern char latex[MAXSTR];
extern char viewdvi[MAXSTR];
extern char Rcommand[MAXSTR];

/* global GUI equipment */
extern windata_t *mdata;
extern GtkTargetEntry gretl_drag_targets[];
extern PangoFontDescription *fixed_font;
#ifdef GRETL_EDIT
extern GtkWidget *editor;
#endif

#include "gretl_enums.h"

/* functions follow */

#ifndef WIN32
void set_wm_icon (GtkWidget *w);
#endif

void set_tryfile (const char *fname);
char *get_tryfile (void);
void clear_tryfile (void);
int tryfile_is_set (void);
int should_ignore_rc (void);
void about_dialog (GtkWidget *w);

#ifdef GRETL_EDIT
void set_editor (GtkWidget *w);
gboolean open_tryfile (gboolean startup);
#else
gboolean open_tryfile (gboolean startup, gboolean dnd);
int mdata_selection_count (void);
int mdata_active_var (void);
void populate_varlist (void);
void clear_varlist (GtkWidget *widget);
void mdata_select_last_var (void);
int gui_restore_sample (DATASET *dset);
void make_list_from_main (void);
void do_stop_script (GtkWidget *w, gpointer p);
int clear_stop_script (PRN *prn);
void show_link_cursor (GtkWidget *w, gpointer p);
gchar *user_friendly_menu_path (const char *mpath,
				gboolean modelwin);
int is_control_key (guint k);
int mainwin_get_vwin_insertion (void);
int mainwin_insert_vwin (windata_t *vwin);
#endif

#endif /* GRETL_H */
