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

#ifndef GRETLTYPES_H
#define GRETLTYPES_H

#define GRETL_STOCK_TEX     "gretl-tex"
#define GRETL_STOCK_MAIL    "gretl-mail"
#define GRETL_STOCK_TS      "gretl-tsplot"
#define GRETL_STOCK_BOX     "gretl-boxplot"
#define GRETL_STOCK_PDF     "gretl-pdf"
#define GRETL_STOCK_BOOK    "gretl-book"
#define GRETL_STOCK_CALC    "gretl-calc"
#define GRETL_STOCK_ICONS   "gretl-icons"
#define GRETL_STOCK_MODEL   "gretl-model"
#define GRETL_STOCK_CONSOLE "gretl-console"
#define GRETL_STOCK_SCATTER "gretl-scatter"
#define GRETL_STOCK_FUNC    "gretl-func"
#define GRETL_STOCK_PIN     "gretl-pin"
#define GRETL_STOCK_ALPHA   "gretl-alpha"
#define GRETL_STOCK_EN      "gretl-en"
#define GRETL_STOCK_SPLIT_H "gretl-split-h"
#define GRETL_STOCK_SPLIT_V "gretl-split-v"


#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 8)
# define NO_INFO_ICON 1
# define GRETL_STOCK_INFO "gretl-info"
# define GTK_STOCK_INFO GRETL_STOCK_INFO
#endif

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 6)
# define NO_EDIT_ICON 1
# define GRETL_STOCK_EDIT "gretl-edit"
# define GTK_STOCK_EDIT GRETL_STOCK_EDIT
# define GRETL_STOCK_SCRIPT "gretl-script"
# define GTK_STOCK_ABOUT NULL
#endif

enum windata_flags {
    VWIN_HELP_ACTIVE     = 1 << 0,
    VWIN_BUSY            = 1 << 1,
    VWIN_DELETE_FNAME    = 1 << 2,
    VWIN_STICKY          = 1 << 3,
    VWIN_CONTENT_CHANGED = 1 << 4,
    VWIN_SESSION_GRAPH   = 1 << 5
};

typedef struct _windata_t windata_t;

#include <gtksourceview/gtksourceview.h>

struct _windata_t {
    GtkWidget *main;      /* top-level GTK window */
    GtkWidget *vbox;      /* vbox within main */
    GtkWidget *text;      /* text or sourceview object */
    GtkWidget *listbox;   /* or: box containing tree or list */
    GtkWidget *mbar;      /* menubar, or toolbar */
    GtkWidget *finder;    /* search entry in top bar */
    GtkWidget *status;    /* status label */
    GtkWidget *popup;     /* popup menu */
    GtkUIManager *ui;     /* UI definition */
    windata_t *gretl_parent;
    windata_t **gretl_children;
    gpointer data;
    int active_var; 
    int role;
    int n_model_tests;
    int n_gretl_children;
    unsigned char flags;
    char fname[MAXLEN];
    GtkSourceBuffer *sbuf;
};

typedef struct dialog_opts_ dialog_opts;

struct dialog_opts_ {
    int n;
    int type;
    gretlopt *optp;
    const gretlopt *vals;
    const char **strs;
};

typedef struct GretlToolItem_ GretlToolItem;

struct GretlToolItem_ {
    const gchar *tip;
    const gchar *icon;
    GCallback func;
    int flag;
};

#define window_is_busy(w)    (w->flags & VWIN_BUSY)

#define window_help_is_active(w)    (w->flags & VWIN_HELP_ACTIVE)
#define set_window_help_active(w)   (w->flags |= VWIN_HELP_ACTIVE)
#define unset_window_help_active(w) (w->flags &= ~VWIN_HELP_ACTIVE)

#define window_delete_filename(w)       (w->flags & VWIN_DELETE_FNAME)
#define set_window_delete_filename(w)   (w->flags |= VWIN_DELETE_FNAME)
#define unset_window_delete_filename(w) (w->flags &= ~VWIN_DELETE_FNAME)

#endif /* GRETLTYPES_H */
