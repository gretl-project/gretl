/*
 *   Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* gretl.h, main header file for gretl gui */

#ifndef GRETL_H
#define GRETL_H

#ifdef _WIN32
# include "winconfig.h"
#else
# include "config.h"
# ifdef USE_GNOME
#   include <gnome.h>
# endif
#endif

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include "libgretl.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#if GTK_MAJOR_VERSION < 2
# define OLD_GTK
# include <gtkextra/gtkiconfilesel.h>
#else
# define GNULL (gconstpointer) NULL
#endif

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
#endif

#include "gretltypes.h"
#include "gui_utils.h"
#include "callbacks.h"
#include "dialogs.h"
#include "library.h"
#include "settings.h"
#include "helpfiles.h"

#ifdef ENABLE_NLS
# include "locale.h"
#endif

#define GRETL_BUFSIZE 8192
#define MAXSTR 255

/* basic global program vars */
extern double **Z;
extern DATAINFO *datainfo;
extern PATHS paths; 
extern MODEL **models;

/* global counters */
extern int plot_count;
extern int orig_vars;

/* global state variables */
extern int data_status;
extern int nls_on;
extern char line[1024];
extern int *default_list;
extern char *storelist;
extern gchar *clipboard_buf;
extern float gui_scale;

/* global filenames */
extern char cmdfile[MAXLEN];
extern char scriptfile[MAXLEN];
extern char trydatfile[MAXLEN];
extern char tryscript[MAXLEN];

/* global error string */
extern char *errtext;

/* global option-related vars */
extern int expert;
extern int updater;
#ifdef G_OS_WIN32
extern int wimp;
#endif
extern char viewdvi[MAXSTR];
extern char viewps[MAXSTR];

/* global GUI equipment */
extern windata_t *mdata;
extern GtkTargetEntry gretl_drag_targets[];
#ifdef OLD_GTK
extern GdkFont *fixed_font;
#else
extern PangoFontDescription *fixed_font;
#endif

#include "gretl_enums.h"

/* functions follow */

#ifndef WIN32
int gretl_fork (const char *prog, const char *arg);
#endif
 
void main_menubar_state (gboolean s);

void edit_info_state (gboolean s);

void populate_varlist (void);

void clear_varlist (GtkWidget *widget);

void clear_sample_label (void);

void set_sample_label (DATAINFO *pdinfo);

void restore_sample_state (gboolean s); 

void restore_sample (void);

void refresh_data (void);

gint main_popup (GtkWidget *widget, GdkEventButton *event, 
		 gpointer data);

/* functions defined in files other than gretl.c */

void file_selector (const char *msg, int action, gpointer data);
void gui_get_series (gpointer data, guint bci_code, GtkWidget *widget);
void import_db_series (windata_t *dbwin);
void display_files (gpointer data, guint code, GtkWidget *widget);
void gpt_save_dialog (void);
void compact_data_set (void);
void stats_calculator (gpointer data, guint code, GtkWidget *widget);

#ifndef WIN32
void set_wm_icon (GtkWidget *w, gpointer data);
#else
int create_child_process (char *prog, char *env);
#endif

#endif /* GRETL_H */
