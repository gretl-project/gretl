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
#define GRETL_STOCK_CP      "gretl-coeffplot"
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
#define GRETL_STOCK_JOIN_H  "gretl-join-h"
#define GRETL_STOCK_JOIN_V  "gretl-join-v"
#define GRETL_STOCK_WINLIST "gretl-winlist"
#define GRETL_STOCK_BUNDLE  "gretl-bundle"
#define GRETL_STOCK_DB      "gretl-db"
#define GRETL_STOCK_GRETL   "gretl-gretl"
#define GRETL_STOCK_TABLE   "gretl-table"
#define GRETL_STOCK_PAGE    "gretl-page"
#define GRETL_STOCK_TOOLS   "gretl-tools"
#define GRETL_STOCK_BIGGER  "gretl-bigger"
#define GRETL_STOCK_SMALLER "gretl-smaller"
#define GRETL_STOCK_MENU    "gretl-menu"
#define GRETL_STOCK_HMAP    "gretl-hmap"
#define GRETL_STOCK_DBN     "gretl-dbnomics"
#define GRETL_STOCK_FCAST   "gretl-fcast"
#define GRETL_STOCK_CLOSE   "gretl-close"
#define GRETL_STOCK_QUERY   "gretl-query"
#define GRETL_STOCK_EYE     "gretl-eye"
#define GRETL_STOCK_EYE_OFF "gretl-eye-off"
#define GRETL_STOCK_MD      "gretl-markdown"

typedef enum {
    VWIN_HELP_ACTIVE     = 1 << 0,
    VWIN_DELETE_FNAME    = 1 << 1,
    VWIN_STICKY          = 1 << 2,
    VWIN_CONTENT_CHANGED = 1 << 3,
    VWIN_SESSION_GRAPH   = 1 << 4,
    VWIN_TABBED          = 1 << 5,
    VWIN_CB_PDF          = 1 << 6,
    VWIN_MULTI_SERIES    = 1 << 7,
    VWIN_NO_SAVE         = 1 << 8,
    VWIN_USE_FOOTER      = 1 << 9,
    WVIN_KEY_SIGNAL_SET  = 1 << 10,
    VWIN_SWALLOW         = 1 << 11
} windata_flags;

typedef struct windata_t_ windata_t;

struct windata_t_ {
    GtkWidget *main;      /* top-level GTK window */
    GtkWidget *topmain;   /* for use when embedded in tabs */
    GtkWidget *hpanes1;   /* upper horizontally opposed panes */
    GtkWidget *hpanes2;   /* lower horizontally opposed panes */
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
    windata_flags flags;
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

#define window_help_is_active(w)    (w->flags & VWIN_HELP_ACTIVE)
#define set_window_help_active(w)   (w->flags |= VWIN_HELP_ACTIVE)
#define unset_window_help_active(w) (w->flags &= ~VWIN_HELP_ACTIVE)
#define window_is_tab(w)            (w->flags & VWIN_TABBED)

#define window_delete_filename(w)     (w->flags & VWIN_DELETE_FNAME)
#define set_window_delete_filename(w) (w->flags |= VWIN_DELETE_FNAME)

#endif /* GRETLTYPES_H */
