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

#ifdef OS_WIN32
# include "../winconfig.h"
#else
# include "../config.h"
# include <gtkextra/gtkiconfilesel.h>
# ifdef USE_GNOME
#   include <gnome.h>
# endif
#endif
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include "../lib/src/libgretl.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "gretltypes.h"
#include "gui_utils.h"
#include "callbacks.h"
#include "dialogs.h"
#include "lib.h"
#include "session.h"

#define MAXSTR 255

extern windata_t *mdata;
extern DATAINFO *datainfo;
extern char *errtext;
extern char cmdfile[MAXLEN];
extern char savefile[MAXLEN];
extern char line[MAXLEN];
extern char scriptfile[MAXLEN];
extern char fontspec[MAXLEN];
extern char expert[6];
extern char updater[6];
extern gchar datalabel[64];
extern PATHS paths; 
extern double *Z;
extern CMD command;
extern print_t cmds;
extern int plot_count;
extern int popup_connected;

extern GtkWidget *sample_label;
extern GdkFont *fixed_font;

extern int oflag, batch, data_file_open;
extern int orig_vars;
extern GENERATE genr; 
extern MODEL **models;
extern SESSION session;
extern session_t rebuild;
extern int *default_list;
extern char *storelist;
extern gchar *clipboard_buf;

enum file_ops {
    OPEN_DATA = 1,
    OPEN_DB,
    OPEN_RATSDB,
    OPEN_SCRIPT,
    OPEN_SAMPLE_SCRIPT,
    OPEN_CSV,
    OPEN_BOX,
    OPEN_SESSION,
    END_OPEN,      /* marker for end of file open section */
    SAVE_DATA,
    SAVE_GZDATA,
    SAVE_BIN1,
    SAVE_BIN2,
    EXPORT_OCTAVE,
    EXPORT_R,
    EXPORT_R_ALT,
    EXPORT_CSV,
    END_SAVE_DATA,  /* marker for end of data-saving section */
    SAVE_CMDS,
    SAVE_TEX_TAB,
    SAVE_TEX_EQ,
    SAVE_SCRIPT,
    SAVE_OUTPUT,
    SAVE_SESSION,
    SAVE_MODEL,
    SAVE_GNUPLOT,
    SAVE_LAST_GRAPH,
    SAVE_GP_CMDS,
    SAVE_CONSOLE,
    SAVE_HTML,
    OP_MAX
};

enum editables {
    EDIT_LOG = 1,
    EDIT_OUTPUT,
    EDIT_SSHEET
};

enum stat_codes {
    ESS = 1,
    R2,
    TR2,
    DF,
    SIGMA
};

enum browser_codes {
    RAMU_DATA = 1,
    RAMU_PS,
    GREENE_DATA,
    GREENE_PS,
    PWT_DATA,
    PWT_PS,
    NATIVE_DB,
    RATS_DB,
    REMOTE_DB,
    NATIVE_SERIES,
    RATS_SERIES,
    REMOTE_SERIES,
    MAINWIN
};

enum exec_codes {
    CONSOLE_EXEC,
    SCRIPT_EXEC,
    SESSION_EXEC
};

enum cgi_options {
    LIST_DBS = 1,
    GRAB_IDX,
    GRAB_DATA,
    SHOW_IDX,
    SHOW_DBS,
    GRAB_NBO_DATA
};
    
enum clipstuff {
    TARGET_STRING,
    TARGET_TEXT,
    TARGET_COMPOUND_TEXT
};

enum copy_variants {
    COPY_SELECTION,
    COPY_TEXT,
    COPY_HTML,
    COPY_LATEX,
    COPY_RTF
};

/* functions follow */
 
void menubar_state (gboolean s);

void graphmenu_state (gboolean s);

void console_state (gboolean s);

gint populate_clist (GtkWidget *widget, DATAINFO *pdatainfo);

void clear_sample_label (void);

void set_sample_label (DATAINFO *pdinfo);

gint main_popup (GtkWidget *widget, GdkEventButton *event, 
		 gpointer data);

#ifndef G_OS_WIN32
void set_wm_icon (GtkWidget *w, gpointer data);
#else
# define isnan(x) ((x) != (x))
#endif

#endif /* GRETL_H */
