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

#include "config.h"
#include <gtkextra/gtkiconfilesel.h>

#ifdef USE_GNOME
# include <gnome.h>
#endif

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "libgretl.h"

#include "gretltypes.h"
#include "gui_utils.h"
#include "callbacks.h"
#include "dialogs.h"
#include "library.h"
#include "session.h"
#include "settings.h"
#include "helpfiles.h"

#define MAXSTR 255

/* basic global program vars */
extern double **Z;
extern DATAINFO *datainfo;
extern PATHS paths; 
extern CMD command;
extern PRN *cmds;
extern MODEL **models;
extern SESSION session;
extern SESSIONBUILD rebuild;

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
extern char viewdvi[MAXSTR];

/* global GUI equipment */
extern windata_t *mdata;
extern GdkFont *fixed_font;
extern GtkTargetEntry gretl_drag_targets[];

enum extra_cmds {
    RENAME = NC,
    RELABEL,
    VSETMISS,
    GSETMISS,
    SMPLDUM,
    SMPLBOOL,
    MARKERS,
    STORE_MODEL,
    VAR_SUMMARY,
    CORR_SELECTED,
    SUMMARY_SELECTED,
    GENR_NORMAL,
    GENR_UNIFORM,
    ONLINE,
    EXPORT,
    MEANTEST2,
    MODEL_GENR,
    GR_PLOT,
    GR_XY,
    GR_IMP,
    GR_DUMMY,
    GR_BOX,
    GR_NBOX,
    COMPACT,
    COEFFINT,
    COVAR,
    STAT_TABLE,
    H_TEST,
    TRAMO,
    X12A,
    RANGE_MEAN,
    VIEW_SERIES,
    VIEW_MODEL,
    VIEW_LOG,
    VIEW_DATA,
    VIEW_SCRIPT,
    VIEW_CODEBOOK,
    DATA_REPORT,
    SCRIPT_OUT,
    CONSOLE,
    EDIT_HEADER,
    EDIT_SCRIPT,
    EDIT_NOTES,
    CLI_HELP,
    GUI_HELP,
    MODELTABLE,
    CMD_MAX
};

enum file_ops {
    OPEN_DATA = CMD_MAX + 1, /* don't collide with extra_cmds */
    OPEN_DB,
    OPEN_RATSDB,
    OPEN_SCRIPT,
    OPEN_CSV,
    APPEND_CSV,
    OPEN_BOX,
    OPEN_GNUMERIC,
    APPEND_GNUMERIC,
    OPEN_EXCEL,
    APPEND_EXCEL,
    OPEN_DES,
    OPEN_SESSION,
    END_OPEN,      /* marker for end of file open section */
    SAVE_DATA,
    SAVE_DATA_AS,
    SAVE_GZDATA,
    SAVE_BIN1,
    SAVE_BIN2,
    EXPORT_OCTAVE,
    EXPORT_R,
    EXPORT_R_ALT,
    EXPORT_CSV,
    COPY_CSV,
    END_SAVE_DATA,  /* marker for end of data-saving section */
    SAVE_CMDS,
    SAVE_TEX_TAB,
    SAVE_TEX_EQ,
    SAVE_TEX_TAB_FRAG,
    SAVE_TEX_EQ_FRAG,
    SAVE_SCRIPT,
    SAVE_OUTPUT,
    SAVE_SESSION,
    SAVE_MODEL,
    SAVE_GNUPLOT,
    SAVE_BOXPLOT_EPS,
    SAVE_BOXPLOT_PS,
    SAVE_BOXPLOT_XPM,
    SAVE_LAST_GRAPH,
    SAVE_THIS_GRAPH,
    SAVE_GP_CMDS,
    SAVE_CONSOLE,
    SET_PATH,
    FILE_OP_MAX
};

#define SAVE_DATA_ACTION(i) (i >= SAVE_DATA && i < END_SAVE_DATA)

enum browser_codes {
    TEXTBOOK_DATA = FILE_OP_MAX + 1, /* don't collide with file_ops enum */
    PS_FILES,
    NATIVE_DB,
    RATS_DB,
    REMOTE_DB,
    NATIVE_SERIES,
    RATS_SERIES,
    REMOTE_SERIES,
    MAINWIN
};

enum stat_codes {
    ESS = 1,
    R2,
    TR2,
    DF,
    SIGMA,
    LNL
};

enum exec_codes {
    CONSOLE_EXEC,
    SCRIPT_EXEC,
    SESSION_EXEC,
    REBUILD_EXEC,
    SAVE_SESSION_EXEC
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
    COPY_SELECTION = 1,
    COPY_TEXT,
    COPY_LATEX,
    COPY_LATEX_EQUATION,
    COPY_RTF
};

enum data_status {
    HAVE_DATA     = 1 << 0,
    BOOK_DATA     = 1 << 1,
    USER_DATA     = 1 << 2,
    IMPORT_DATA   = 1 << 3,
    GUI_DATA      = 1 << 4,
    MODIFIED_DATA = 1 << 5
};

enum drag_types {
    GRETL_FILENAME,
    GRETL_POINTER,
    GRETL_MODEL_POINTER
};

enum file_lists {
    FILE_LIST_DATA = 1,
    FILE_LIST_SESSION,
    FILE_LIST_SCRIPT,
};

enum graph_types {
    GRETL_GNUPLOT_GRAPH,
    GRETL_BOXPLOT
};

/* functions follow */

void gretl_fork (const char *prog, const char *arg);
 
void menubar_state (gboolean s);

#ifndef GNUPLOT_PNG
void graphmenu_state (gboolean s);
#endif

gint populate_main_varlist (void);

void clear_sample_label (void);

void set_sample_label (DATAINFO *pdinfo);

void restore_sample_state (gboolean s); 

void restore_sample (gpointer data, int verbose, GtkWidget *w);

void refresh_data (void);

gint main_popup (GtkWidget *widget, GdkEventButton *event, 
		 gpointer data);

/* functions defined in files other than gretl.c */

void file_selector (const char *msg, int action, gpointer data);
void gui_get_series (gpointer data, guint bci_code, 
		     GtkWidget *widget);
void import_db_series (windata_t *dbwin);
void display_files (gpointer data, guint code, GtkWidget *widget);
void gpt_save_dialog (void);
void compact_data_set (void);

/* webget.c */
int update_query (int verbose); 
int retrieve_url (int opt, const char *dbase, const char *series, 
		  int filesave, char **saver, char *errbuf);
int proxy_init (const char *dbproxy);

void set_wm_icon (GtkWidget *w, gpointer data);

#endif /* GRETL_H */
