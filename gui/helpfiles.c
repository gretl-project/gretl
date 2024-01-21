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

/* helpfiles.c for gretl */

#include "gretl.h"
#include "textbuf.h"
#include "gretl_www.h"
#include "dlgutils.h"
#include "toolbar.h"
#include "winstack.h"
#include "base_utils.h"

#ifndef GRETL_EDIT
#include "treeutils.h"
#include "menustate.h"
#include "database.h"
#include "fncall.h"
#endif

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
# include <dirent.h>
#endif

#ifdef OS_OSX
# include "osx_open.h"
#endif

#define HDEBUG 0

#define cmdref_role(r)  (r == CMD_HELP || r == CMD_HELP_EN)
#define funcref_role(r) (r == FUNC_HELP || r == FUNC_HELP_EN)

static int translated_cmdref;
static int translated_fnref;

static windata_t *real_do_help (int idx, int pos, int role);
static void en_help_callback (GtkWidget *w, windata_t *hwin);
static void helpwin_set_topic_index (windata_t *hwin, int idx);

/* searching stuff */
static void find_in_text (GtkWidget *button, GtkWidget *dialog);
static void find_string_dialog (void (*findfunc)(), windata_t *vwin);
static gboolean real_find_in_text (GtkTextView *view, const gchar *s,
				   gboolean sensitive,
				   gboolean from_cursor,
				   gboolean search_all);
static void vwin_finder_callback (GtkEntry *entry, windata_t *vwin);
#ifndef GRETL_EDIT
static void find_in_listbox (GtkWidget *button, GtkWidget *dialog);
static gboolean real_find_in_listbox (windata_t *vwin, const gchar *s,
				      gboolean sensitive,
				      gboolean vnames);
#endif

static GtkWidget *find_dialog = NULL;
static GtkWidget *find_entry;
static char *needle;

enum {
    /* don't collide with the enumeration in toolbar.c */
    EN_ITEM = 255
};

static GretlToolItem help_tools[] = {
    { N_("Larger"), GTK_STOCK_ZOOM_IN, G_CALLBACK(text_larger), 0},
    { N_("Smaller"), GTK_STOCK_ZOOM_OUT, G_CALLBACK(text_smaller), 0},
    { N_("Show English help"), GRETL_STOCK_EN, G_CALLBACK(en_help_callback), EN_ITEM }
};

static gint n_help_tools = G_N_ELEMENTS(help_tools);

struct gui_help_item {
    int code;
    char *string;
};

/* codes and strings for GUI help items other than
   regular gretl commands */

static struct gui_help_item gui_help_items[] = {
    { 0,              "nothing" },
    { GR_PLOT,        "graphing" },
    { GR_XY,          "graphing" },
    { GR_DUMMY,       "factorized" },
    { GR_XYZ,         "controlled" },
    { BXPLOT,         "boxplots" },
    { GR_BOX,         "boxplots" },
    { GR_FBOX,        "boxplots" },
    { GR_3D,          "3-D" },
    { ONLINE,         "online" },
    { EXPORT,         "export" },
    { SMPLBOOL,       "sampling" },
    { SMPLDUM,        "sampling" },
    { COMPACT,        "compact" },
    { EXPAND,         "expand" },
    { TDISAGG,        "tdisagg" },
    { VSETMISS,       "missing" },
    { GSETMISS,       "missing" },
    { GUI_HELP,       "dialog" },
    { GENR_RANDOM,    "genrand" },
    { SEED_RANDOM,    "genseed" },
    { KERNEL_DENSITY, "density" },
    { HCCME,          "hccme" },
    { IRF_BOOT,       "irfboot" },
    { HTEST,          "gui-htest" },
    { HTESTNP,        "gui-htest-np" },
    { MODEL_RESTR,    "restrict-model" },
    { SYS_RESTR,      "restrict-system" },
    { VECM_RESTR,     "restrict-vecm" },
    { LAGS_DIALOG,    "lags-dialog" },
    { MINIBUF,        "minibuffer" },
    { ALAGSEL,        "ARMA-lagselect" },
    { VLAGSEL,        "VAR-lagselect" },
    { VAROMIT,        "VAR-omit" },
    { PANEL_MODE,     "panel-mode" },
    { TS_TO_PANEL,    "ts-to-panel" },
    { PANEL_WLS,      "panel-wls" },
    { PANEL_B,        "panel-between" },
    { BOOTSTRAP,      "bootstrap" },
    { TRANSPOS,       "transpos" },
    { DATASORT,       "datasort" },
    { WORKDIR,        "workdir" },
    { DFGLS,          "dfgls" },
    { GPT_ADDLINE,    "addline" },
    { GPT_CURVE,      "curve" },
    { SAVE_SESSION,   "save-session" },
    { SAVE_CMD_LOG,   "save-script" },
    { BFGS_CONFIG,    "bfgs-config" },
    { SAVE_LABELS,    "save-labels" }, /* FIXME */
    { OPEN_LABELS,    "add-labels" },  /* FIXME */
    { COUNTMOD,       "count-model" },
    { REGLS,          "regls" },
    { REGLS_ADV,      "regls-advanced" },
    { BWFILTER,       "bwfilter" },
    { POLYWEIGHTS,    "polyweights" },
    { EMAFILTER,      "ema-filter" },
    { X12AHELP,       "x12a" },
    { MAILHELP,       "mailer" },
    { LOESS,          "loess" },
    { NADARWAT,       "nadarwat" },
    { CLUSTER,        "cluster" },
    { GUI_FUNCS,      "gui-funcs" },
    { MENU_ATTACH,    "menu-attach" },
    { REPROBIT,       "reprobit" },
    { DAILY_PURGE,    "daily-purge" },
    { PKG_FILES,      "data-files" },
    { PKG_DEPS,       "pkg-depends" },
    { EDITOR,         "script-editor" },
    { MIDAS_LIST,     "MIDAS_list" },
    { MIDAS_PARM,     "MIDAS_parm" },
    { PKGHELP,        "packages" },
    { DBNHELP,        "dbnomics" },
    { MAPHELP,        "maps" },
    { MDHELP,         "mdhelp" },
    { KALMAN,         "kalman" },
    { EDITOR,         "gretl_edit" },
    { -1,             NULL },
};

/* state the topic headings from the script help files so they
   can be translated */

const char *intl_topics[] = {
    N_("Dataset"),
    N_("Estimation"),
    N_("Graphs"),
    N_("Prediction"),
    N_("Printing"),
    N_("Programming"),
    N_("Statistics"),
    N_("Tests"),
    N_("Transformations"),
    N_("Utilities")
};

/* Handle non-uniqueness of map from 'extra' command words to codes
   (e.g. both GR_PLOT and GR_XY correspond to "graphing").  We want
   the first code, to find the right place in the help file
   when responding to a "context help" request.
*/

static int gui_ci_to_index (int ci)
{
    if (ci < NC) {
	/* regular gretl command, no problem */
	return ci;
    } else {
	int i, k, ret = ci;

	for (i=1; gui_help_items[i].code > 0; i++) {
	    if (ci == gui_help_items[i].code) {
		for (k=i-1; k>0; k--) {
		    /* back up the list so long as the word above
		       is the same as the current one */
		    if (!strcmp(gui_help_items[k].string,
				gui_help_items[i].string)) {
			ret = gui_help_items[k].code;
		    } else {
			break;
		    }
		}
		return ret;
	    }
	}
    }

    return -1;
}

/* Public because it's wanted in textbuf.c for getting
   the help indices of certain items that are not actual
   gretl commands.
*/

int extra_command_number (const char *s)
{
    int i;

    for (i=1; gui_help_items[i].code > 0; i++) {
	if (!strcmp(s, gui_help_items[i].string)) {
	    return gui_help_items[i].code;
	}
    }

    return -1;
}

void helpfile_init (void)
{
    translated_cmdref = using_translated_helpfile(GRETL_CMDREF);
    translated_fnref = using_translated_helpfile(GRETL_FUNCREF);
}

char *quoted_help_string (const char *s)
{
    const char *p, *q;

    p = strchr(s, '"');
    q = strrchr(s, '"');

    if (p != NULL && q != NULL && q - p > 1) {
	p++;
	return g_strndup(p, q - p);
    }

    return g_strdup("Missing string");
}

static int gui_help_topic_index (const char *word)
{
    int h = gretl_command_number(word);

    if (h == 0) {
	h = extra_command_number(word);
    }

    return h;
}

enum {
    STRING_COL,
    POSITION_COL,
    INDEX_COL,
    NUM_COLS
};

static void help_tree_select_row (GtkTreeSelection *selection,
				  windata_t *hwin)
{
    GtkTreeIter iter;
    GtkTreeModel *model;

    if (!gtk_tree_selection_get_selected(selection, &model, &iter)) {
	return;
    }

    if (!gtk_tree_model_iter_has_child(model, &iter)) {
	int pos, idx;

	gtk_tree_model_get(model, &iter,
			   POSITION_COL, &pos,
			   INDEX_COL, &idx,
			   -1);

	if (idx != hwin->active_var) {
	    /* not already in position */
	    real_do_help(idx, pos, hwin->role);
	}
    }
}

static int get_section_iter (GtkTreeModel *model, const char *s,
			     GtkTreeIter *iter)
{
    gchar *sect;
    int found = 0;

    if (!gtk_tree_model_get_iter_first(model, iter)) {
	return 0;
    }

    while (!found && gtk_tree_model_iter_next(model, iter)) {
	if (gtk_tree_model_iter_has_child(model, iter)) {
	    gtk_tree_model_get(model, iter, STRING_COL, &sect, -1);
	    if (!strcmp(s, sect)) {
		found = 1;
	    }
	    g_free(sect);
	}
    }

    return found;
}

static const char *real_funcs_heading (const char *s)
{
    if (!strcmp(s, "access")) {
	return _("Accessors");
    } else if (!strcmp(s, "straccess")) {
	return _("Built-in strings");
    } else if (!strcmp(s, "data-utils")) {
	return _("Data utilities");
    } else if (!strcmp(s, "math")) {
	return _("Mathematical");
    } else if (!strcmp(s, "transforms")) {
	return _("Transformations");
    } else if (!strcmp(s, "matrix")) {
	return _("Matrix manipulation");
    } else if (!strcmp(s, "linalg")) {
	return _("Linear algebra");
    } else if (!strcmp(s, "complex")) {
	return _("Complex numbers");
    } else if (!strcmp(s, "numerical")) {
	return _("Numerical methods");
    } else if (!strcmp(s, "probdist")) {
	return _("Probability");
    } else if (!strcmp(s, "panel")) {
	return _("Panel data");
    } else if (!strcmp(s, "calendar")) {
	return _("Calendar");
    } else if (!strcmp(s, "timeseries")) {
	return _("Time-series");
    } else if (!strcmp(s, "stats")) {
	return _("Statistical");
    } else if (!strcmp(s, "nonparam")) {
	return _("Non-parametric");
    } else if (!strcmp(s, "midas")) {
	return _("MIDAS");
    } else if (!strcmp(s, "sspace")) {
	return _("State space");
    } else if (!strcmp(s, "programming")) {
	return _("Programming");
    } else if (!strcmp(s, "strings")) {
	return _("Strings");
    } else if (!strcmp(s, "mpi")) {
	return _("MPI");
    } else {
	return "??";
    }
}

static GtkTreeStore *make_help_topics_tree (int role)
{
    const char *fname;
    GtkTreeStore *store;
    GtkTreeIter iter, parent;
    gchar *s, *buf = NULL;
    char word[32], sect[48];
    const char *heading;
    int pos = 0, idx = 0;
    int err;

    if (role == CMD_HELP) {
	fname = helpfile_path(GRETL_CMDREF, 0, 0);
    } else if (role == GUI_HELP) {
	fname = helpfile_path(GRETL_GUI_HELP, 0, 0);
    } else if (role == FUNC_HELP) {
	fname = helpfile_path(GRETL_FUNCREF, 0, 0);
    } else if (role == CMD_HELP_EN) {
	fname = helpfile_path(GRETL_CMDREF, 0, 1);
    } else if (role == GUI_HELP_EN) {
	fname = helpfile_path(GRETL_GUI_HELP, 0, 1);
    } else if (role == FUNC_HELP_EN) {
	fname = helpfile_path(GRETL_FUNCREF, 0, 1);
    } else {
	return NULL;
    }

    err = gretl_file_get_contents(fname, &buf, NULL);
    if (err) {
	return NULL;
    }

    store = gtk_tree_store_new(NUM_COLS, G_TYPE_STRING,
			       G_TYPE_INT, G_TYPE_INT);

    gtk_tree_store_append(store, &iter, NULL);
    gtk_tree_store_set(store, &iter, STRING_COL, _("Index"),
		       POSITION_COL, 0, -1);

    s = buf;

#if HDEBUG
    fprintf(stderr, "*** processing %s ***\n", fname);
#endif

    while (*s) {
	if (*s == '\n' && *(s+1) == '#' &&
	    *(s+2) != '#' && *(s+2) != '\0') {
	    if (sscanf(s+2, "%31s %47s", word, sect) == 2) {
		if (funcref_role(role)) {
		    heading = real_funcs_heading(sect);
		} else {
		    heading = _(sect);
		}
		if (!get_section_iter(GTK_TREE_MODEL(store), heading, &parent)) {
		    gtk_tree_store_append(store, &parent, NULL);
		    gtk_tree_store_set(store, &parent,
				       STRING_COL, heading,
				       POSITION_COL, 0,
				       INDEX_COL, 0,
				       -1);
		}
		gtk_tree_store_append(store, &iter, &parent);
		if (funcref_role(role)) {
		    ++idx;
		} else if (role == GUI_HELP || role == GUI_HELP_EN) {
		    idx = gui_help_topic_index(word);
		} else {
		    idx = gretl_help_index(word);
		}
#if HDEBUG
		fprintf(stderr, " %s: pos %d, idx %d\n", word, pos+1, idx);
#endif
		gtk_tree_store_set(store, &iter,
				   STRING_COL, word,
				   POSITION_COL, pos + 1,
				   INDEX_COL, idx,
				   -1);
	    }
	}
	s++;
	pos++;
    }

    g_free(buf);

    return store;
}

static GtkTreeStore *get_help_topics_tree (int role)
{
    static GtkTreeStore *cmd_tree;
    static GtkTreeStore *gui_tree;
    static GtkTreeStore *func_tree;
    static GtkTreeStore *en_cmd_tree;
    static GtkTreeStore *en_gui_tree;
    static GtkTreeStore *en_func_tree;
    GtkTreeStore **ptree = NULL;

    if (role == CMD_HELP) {
	ptree = &cmd_tree;
    } else if (role == GUI_HELP) {
	ptree = &gui_tree;
    } else if (role == FUNC_HELP) {
	ptree = &func_tree;
    } else if (role == CMD_HELP_EN) {
	ptree = &en_cmd_tree;
    } else if (role == GUI_HELP_EN) {
	ptree = &en_gui_tree;
    } else if (role == FUNC_HELP_EN) {
	ptree = &en_func_tree;
    } else {
	return NULL;
    }

    if (*ptree == NULL) {
	*ptree = make_help_topics_tree(role);
    }

    return *ptree;
}

/* add a tree-style navigation pane to the left of the help index or
   text */

int add_help_navigator (windata_t *vwin, GtkWidget *hp)
{
    GtkTreeStore *store;
    GtkWidget *view, *sw;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;

    store = get_help_topics_tree(vwin->role);
    if (store == NULL) {
	return 1;
    }

    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);
    g_object_set(view, "enable-tree-lines", TRUE, NULL);

    renderer = gtk_cell_renderer_text_new();

    column = gtk_tree_view_column_new_with_attributes("",
						      renderer,
						      "text", 0,
						      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
    g_signal_connect(G_OBJECT(select), "changed",
		     G_CALLBACK(help_tree_select_row),
		     vwin);

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(sw), view);
    gtk_paned_pack1(GTK_PANED(hp), sw, FALSE, TRUE);
    gtk_widget_set_size_request(sw, 150, -1);
    gtk_tree_view_columns_autosize(GTK_TREE_VIEW(view));

    g_object_set_data(G_OBJECT(vwin->text), "tview", view);

    return 0;
}

static int help_attr_from_word (const char *word,
				int role, int col)
{
    GtkTreeModel *model;
    GtkTreeIter iter, child;
    gchar *s;
    int attr = 0;

    model = GTK_TREE_MODEL(get_help_topics_tree(role));

    if (!model || !gtk_tree_model_get_iter_first(model, &iter)) {
	return 0;
    }

    while (gtk_tree_model_iter_next(model, &iter)) {
	if (gtk_tree_model_iter_children(model, &child, &iter)) {
	    while (1) {
		gtk_tree_model_get(model, &child, STRING_COL, &s, -1);
		if (!strcmp(s, word)) {
		    gtk_tree_model_get(model, &child, col, &attr, -1);
		    g_free(s);
		    return attr;
		}
		g_free(s);
		if (!gtk_tree_model_iter_next(model, &child)) {
		    break;
		}
	    }
	}
    }

    return 0;
}

int function_help_index_from_word (const char *word, int role)
{
    return help_attr_from_word(word, role, INDEX_COL);
}

static int function_help_pos_from_word (const char *word, int role)
{
    return help_attr_from_word(word, role, POSITION_COL);
}

static int help_pos_from_index (int idx, int role)
{
    GtkTreeModel *model;
    GtkTreeIter iter, child;
    int pos, tidx;

    model = GTK_TREE_MODEL(get_help_topics_tree(role));

    if (!model || !gtk_tree_model_get_iter_first(model, &iter)) {
	return 0;
    }

    while (gtk_tree_model_iter_next(model, &iter)) {
	if (gtk_tree_model_iter_children(model, &child, &iter)) {
	    while (1) {
		gtk_tree_model_get(model, &child, INDEX_COL, &tidx, -1);
		if (tidx == idx) {
		    gtk_tree_model_get(model, &child, POSITION_COL, &pos, -1);
		    return pos;
		}
		if (!gtk_tree_model_iter_next(model, &child)) {
		    break;
		}
	    }
	}
    }

    return 0;
}

static int help_index_from_pos (int pos, int role)
{
    GtkTreeModel *model;
    GtkTreeIter iter, child;
    int idx, tpos;

    model = GTK_TREE_MODEL(get_help_topics_tree(role));

    if (!model || !gtk_tree_model_get_iter_first(model, &iter)) {
	return 0;
    }

    while (gtk_tree_model_iter_next(model, &iter)) {
	if (gtk_tree_model_iter_children(model, &child, &iter)) {
	    while (1) {
		gtk_tree_model_get(model, &child, POSITION_COL, &tpos, -1);
		if (tpos == pos) {
		    gtk_tree_model_get(model, &child, INDEX_COL, &idx, -1);
		    return idx;
		}
		if (!gtk_tree_model_iter_next(model, &child)) {
		    break;
		}
	    }
	}
    }

    return 0;
}

static gboolean help_iter_from_index (int idx, int role,
				      GtkTreeIter *iter,
				      GtkTreeIter *parent)
{
    GtkTreeModel *model;
    int fnum;

    model = GTK_TREE_MODEL(get_help_topics_tree(role));

    if (!model || !gtk_tree_model_get_iter_first(model, parent)) {
	return 0;
    }

    while (gtk_tree_model_iter_next(model, parent)) {
	if (gtk_tree_model_iter_children(model, iter, parent)) {
	    while (1) {
		gtk_tree_model_get(model, iter, INDEX_COL, &fnum, -1);
		if (idx == fnum) {
		    return TRUE;
		}
		if (!gtk_tree_model_iter_next(model, iter)) {
		    break;
		}
	    }
	}
    }

    return FALSE;
}

static void en_help_callback (GtkWidget *w, windata_t *hwin)
{
    int idx = hwin->active_var;
    int pos, role = 0;

    if (hwin->role == CMD_HELP) {
	role = CMD_HELP_EN;
    } else if (hwin->role == GUI_HELP) {
	role = GUI_HELP_EN;
    } else if (hwin->role == FUNC_HELP) {
	role = FUNC_HELP_EN;
    } else {
	return;
    }

    pos = help_pos_from_index(idx, role);

    if (pos < 0 && role != GUI_HELP_EN) {
	/* missed, but we can at least show the index page */
	pos = 0;
    }

    real_do_help(idx, pos, role);
}

#if GTK_MAJOR_VERSION == 2

static void normalize_base (GtkWidget *w, gpointer p)
{
    gtk_widget_modify_base(w, GTK_STATE_SELECTED, NULL);
}

void notify_string_not_found (GtkWidget *entry)
{
    GdkColor color;

    gdk_color_parse("red", &color);
    gtk_widget_grab_focus(entry);
    gtk_editable_select_region(GTK_EDITABLE(entry), 0, -1);
    gtk_widget_modify_base(entry, GTK_STATE_SELECTED, &color);
    g_signal_connect(G_OBJECT(entry), "changed",
		     G_CALLBACK(normalize_base), NULL);
}

#else /* GTK 3.0 */

static void normalize_style (GtkWidget *w, gpointer p)
{
    GtkStyleContext *context;

    if (p == NULL) {
	context = gtk_widget_get_style_context(w);
    } else {
	context = GTK_STYLE_CONTEXT(p);
    }
    gtk_style_context_remove_class(context, GTK_STYLE_CLASS_ERROR);
}

void notify_string_not_found (GtkWidget *entry)
{
    GtkStyleContext *context;

    context = gtk_widget_get_style_context(entry);
    gtk_style_context_add_class(context,
				GTK_STYLE_CLASS_ERROR);
    gtk_widget_grab_focus(entry);
    gtk_editable_select_region(GTK_EDITABLE(entry), 0, -1);
    g_signal_connect(G_OBJECT(entry), "changed",
		     G_CALLBACK(normalize_style), context);
}

#endif

#define help_index_ok(r) (r == CMD_HELP || \
                          r == CMD_HELP_EN || \
                          r == FUNC_HELP || \
			  r == FUNC_HELP_EN)

static gboolean finder_key_handler (GtkEntry *entry, GdkEventKey *event,
				    windata_t *vwin)
{
    guint keyval = event->keyval;

#ifdef OS_OSX
    if (keyval == GDK_g && cmd_key(event)) {
	/* Command-G: repeat search */
	vwin_finder_callback(entry, vwin);
	return TRUE;
    }
#endif

    if (keyval == GDK_Tab && help_index_ok(vwin->role) &&
	vwin->active_var == 0) {
	/* tab-completion in help index mode */
	const gchar *s = gtk_entry_get_text(entry);

	if (s != NULL && *s != '\0') {
	    const char *comp = NULL;

	    if (cmdref_role(vwin->role)) {
		comp = gretl_command_complete(s);
	    } else if (funcref_role(vwin->role)) {
		comp = gretl_function_complete(s);
	    }

	    if (comp != NULL) {
		gtk_entry_set_text(entry, comp);
		gtk_editable_set_position(GTK_EDITABLE(entry), -1);
	    }

	    return TRUE;
	}
    } else if (keyval == GDK_g && (event->state & GDK_CONTROL_MASK)) {
	/* Ctrl-G: repeat search */
	vwin_finder_callback(entry, vwin);
	return TRUE;
    }

    return FALSE;
}

#define starts_topic(s) (s[0]=='\n' && s[1]=='#' && s[2]==' ')

/* apparatus for permitting the user to search across all the
   "pages" in a help file
*/

static int maybe_switch_page (const char *s, windata_t *hwin)
{
    const gchar *src, *hbuf;
    int currpos, newpos = 0;
    int wrapped = 0;
    int k, n = strlen(s);
    int ok = 0;

    /* where are we in the help file right now? */
    currpos = help_pos_from_index(hwin->active_var, hwin->role);
    hbuf = (const gchar *) hwin->data;
    k = currpos;
    src = hbuf + k;

    /* skip to start of next page */
    while (*src != '\0') {
	if (starts_topic(src)) {
	    break;
	}
	k++;
	src++;
    }

 retry:

    /* see if the search text can be found on a page other than
       the current one; if so, we'll switch to it */

    while (!ok && *src != '\0') {
	if (starts_topic(src)) {
	    /* record page position */
	    newpos = k + 1;
	} else if (wrapped && newpos == currpos) {
	    /* we got back to where we started */
	    break;
	} else if (newpos != currpos && !strncmp(s, src, n)) {
	    /* found the text on a different page */
	    ok = 1;
	}
	k++;
	src++;
    }

    if (!ok && !wrapped) {
	/* start from the top */
	src = hbuf;
	newpos = k = 0;
	wrapped = 1;
	goto retry;
    }

    if (ok) {
	/* text found: move to new position */
	int idx = help_index_from_pos(newpos, hwin->role);

	set_help_topic_buffer(hwin, newpos);
	helpwin_set_topic_index(hwin, idx);
    }

    return ok;
}

static int all_lower_case (const char *s)
{
    while (*s) {
	if (tolower(*s) != *s) {
	    return 0;
	}
	s++;
    }

    return 1;
}

static int maybe_go_to_page (windata_t *vwin)
{
    gboolean ret = FALSE;
    int idx, pos = 0;

    if (funcref_role(vwin->role)) {
	/* looking for a function */
	idx = function_help_index_from_word(needle, vwin->role);
	if (idx > 0) {
	    pos = help_pos_from_index(idx, vwin->role);
	}
    } else {
	/* looking for a command */
	idx = gretl_command_number(needle);
	if (idx > 0) {
	    pos = help_pos_from_index(idx, vwin->role);
	}
    }

    if (pos > 0) {
	real_do_help(idx, pos, vwin->role);
	ret = TRUE;
    }

    return ret;
}

/* respond to Enter key in the 'finder' entry */

static void vwin_finder_callback (GtkEntry *entry, windata_t *vwin)
{
    gboolean search_all = FALSE;
    gboolean found = FALSE;
    gboolean sensitive;

    needle = gtk_editable_get_chars(GTK_EDITABLE(entry), 0, -1);
    if (needle == NULL || *needle == '\0') {
	return;
    }

    sensitive = !all_lower_case(needle);

#ifdef GRETL_EDIT
    if (g_object_get_data(G_OBJECT(entry), "search-all")) {
	search_all = TRUE;
    }
#else
    if (g_object_get_data(G_OBJECT(entry), "search-all")) {
	if (vwin->role == DBNOMICS_TOP ||
	    vwin->role == DBNOMICS_DB ||
	    vwin->role == DBNOMICS_SERIES) {
	    dbnomics_search(needle, vwin);
	    return;
	} else {
	    search_all = TRUE;
	}
    }
#endif

    if (vwin->text != NULL) {
	gboolean from_cursor = TRUE;

	found = real_find_in_text(GTK_TEXT_VIEW(vwin->text), needle,
				  sensitive, from_cursor, search_all);
    }
#ifndef GRETL_EDIT
    else {
	found = real_find_in_listbox(vwin, needle, sensitive, 0);
    }
#endif

    if (!found && (vwin->role == CMD_HELP ||
		   vwin->role == CMD_HELP_EN ||
		   vwin->role == FUNC_HELP ||
		   vwin->role == FUNC_HELP_EN)) {
	found = maybe_go_to_page(vwin);
    }

    if (!found && search_all) {
	if (maybe_switch_page(needle, vwin)) {
	    found = real_find_in_text(GTK_TEXT_VIEW(vwin->text), needle,
				      sensitive, TRUE, TRUE);
	}
    }

    if (!found) {
	notify_string_not_found(GTK_WIDGET(entry));
    }
}

static void toggle_search_all_help (GtkComboBox *box, GtkWidget *entry)
{
    gint i = gtk_combo_box_get_active(box);

    if (i > 0) {
	g_object_set_data(G_OBJECT(entry), "search-all", GINT_TO_POINTER(1));
    } else {
	g_object_steal_data(G_OBJECT(entry), "search-all");
    }
}

static void finder_add_options (GtkWidget *hbox, GtkWidget *entry)
{
    GtkWidget *combo = gtk_combo_box_text_new();

    combo_box_append_text(combo, _("this page"));
    combo_box_append_text(combo, _("all pages"));
    gtk_box_pack_end(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);

    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(toggle_search_all_help),
		     entry);
}

static void toggle_search_this_help (GtkComboBox *box, GtkWidget *entry)
{
    gint i = gtk_combo_box_get_active(box);

    if (i > 0) {
	g_object_steal_data(G_OBJECT(entry), "search-all");
    } else {
	g_object_set_data(G_OBJECT(entry), "search-all", GINT_TO_POINTER(1));
    }
}

static void finder_add_dbn_options (windata_t *vwin,
				    GtkWidget *hbox,
				    GtkWidget *entry)
{
    GtkWidget *combo = gtk_combo_box_text_new();

    if (vwin->role == DBNOMICS_DB) {
	combo_box_append_text(combo, _("selected dataset"));
    } else if (vwin->role == DBNOMICS_SERIES) {
	combo_box_append_text(combo, _("this dataset"));
    } else {
	combo_box_append_text(combo, _("all DB.NOMICS"));
    }
    combo_box_append_text(combo, _("this window"));
    gtk_box_pack_end(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);

    g_object_set_data(G_OBJECT(entry), "search-all", GINT_TO_POINTER(1));
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(toggle_search_this_help),
		     entry);
}

static void finder_icon_press (GtkEntry *entry,
			       GtkEntryIconPosition pos,
			       GdkEvent *event,
			       windata_t *vwin)
{
    vwin_finder_callback(entry, vwin);
}

static void add_finder_icon (windata_t *vwin, GtkWidget *entry)
{
    gtk_entry_set_icon_from_stock(GTK_ENTRY(entry),
				  GTK_ENTRY_ICON_SECONDARY,
				  GTK_STOCK_FIND);
    gtk_entry_set_icon_activatable(GTK_ENTRY(entry),
				   GTK_ENTRY_ICON_SECONDARY,
				   TRUE);
    g_signal_connect(G_OBJECT(entry), "icon-press",
		     G_CALLBACK(finder_icon_press),
		     vwin);
}

/* add a "search box" to the right of a viewer window's toolbar */

void vwin_add_finder (windata_t *vwin)
{
    GtkWidget *entry;
    GtkWidget *hbox;
    int fwidth = 16;

    hbox = gtk_widget_get_parent(vwin->mbar);
    vwin->finder = entry = gtk_entry_new();

    if (vwin->role == CMD_HELP ||
	vwin->role == CMD_HELP_EN ||
	vwin->role == FUNC_HELP ||
	vwin->role == FUNC_HELP_EN) {
	finder_add_options(hbox, entry);
    } else if (vwin->role == DBNOMICS_TOP ||
	       vwin->role == DBNOMICS_DB ||
	       vwin->role == DBNOMICS_SERIES) {
	finder_add_dbn_options(vwin, hbox, entry);
	fwidth = 28;
    }

    gtk_entry_set_width_chars(GTK_ENTRY(entry), fwidth);
    gtk_box_pack_end(GTK_BOX(hbox), entry, FALSE, FALSE, 5);
    add_finder_icon(vwin, entry);

#ifndef GRETL_EDIT
    if (vwin->role == DBNOMICS_TOP ||
	vwin->role == VIEW_DBSEARCH ||
	vwin->role == DBNOMICS_SERIES ||
	vwin->role == DBNOMICS_DB) {
	maybe_fill_dbn_finder(vwin->finder);
    }
#endif

    g_signal_connect(G_OBJECT(entry), "key-press-event",
		     G_CALLBACK(finder_key_handler), vwin);
    g_signal_connect(G_OBJECT(entry), "activate",
		     G_CALLBACK(vwin_finder_callback),
		     vwin);
}

static void add_footer_close_button (GtkWidget *hbox)
{
    GtkWidget *img = gtk_image_new_from_stock(GRETL_STOCK_CLOSE,
					      GTK_ICON_SIZE_MENU);
    GtkWidget *button = gtk_button_new();

    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    gtk_container_add(GTK_CONTAINER(button), img);
    gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect_swapped(button, "clicked",
			     G_CALLBACK(gtk_widget_hide),
			     hbox);
    gtk_widget_show_all(button);
}

static gint catch_footer_key (GtkWidget *w, GdkEventKey *event,
			      GtkWidget *targ)
{
    if (event->keyval == GDK_Escape) {
	gtk_widget_hide(targ);
	return TRUE;
    } else {
	return FALSE;
    }
}

/* Callback from "hide" signal on the footer finder: we
   want to turn the focus back onto the associated
   text widget
*/

static void vwin_refocus_text (GtkWidget *w, windata_t *vwin)
{
    /* in case the prior string was not found, cancel the
       error indicator
    */
#if GTK_MAJOR_VERSION == 2
    normalize_base(vwin->finder, NULL);
#else
    normalize_style(vwin->finder, NULL);
#endif

    if (vwin_is_editing(vwin)) {
	/* let the text area regain focus */
	gtk_widget_grab_focus(vwin->text);
	gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->text), TRUE);
    }
}

static void vwin_add_footer_finder (windata_t *vwin)
{
    GtkWidget *hbox, *entry;

    hbox = gtk_hbox_new(FALSE, 5);
    vwin->finder = entry = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 20);
    add_finder_icon(vwin, entry);

    gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 10);
    add_footer_close_button(hbox);
    gtk_box_pack_end(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 2);

    g_signal_connect(G_OBJECT(entry), "key-press-event",
		     G_CALLBACK(catch_footer_key), hbox);
    g_signal_connect(G_OBJECT(entry), "key-press-event",
		     G_CALLBACK(finder_key_handler), vwin);
    g_signal_connect(G_OBJECT(entry), "activate",
		     G_CALLBACK(vwin_finder_callback),
		     vwin);
    g_signal_connect(G_OBJECT(hbox), "hide",
		     G_CALLBACK(vwin_refocus_text),
		     vwin);

    gtk_widget_show_all(hbox);
    gtk_widget_grab_focus(entry);
}

#define SHOW_FINDER(r)    (r != GUI_HELP && r != GUI_HELP_EN)
#define SHOW_EN_BUTTON(r) ((translated_cmdref && (r==CMD_HELP || r==GUI_HELP)) || \
			   (translated_fnref && r == FUNC_HELP))

void set_up_helpview_menu (windata_t *hwin)
{
    GretlToolItem *item = NULL;
    GtkWidget *hbox;
    int i;

    hbox = gtk_hbox_new(FALSE, 0);
    hwin->mbar = gretl_toolbar_new(NULL);

    for (i=0; i<n_help_tools; i++) {
	item = &help_tools[i];
	if (!SHOW_EN_BUTTON(hwin->role) && item->flag == EN_ITEM) {
	    continue;
	}
	gretl_toolbar_insert(hwin->mbar, item, item->func, hwin, -1);
    }

    gtk_box_pack_start(GTK_BOX(hwin->vbox), hbox, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), hwin->mbar, FALSE, FALSE, 0);

    vwin_add_winlist(hwin);

    if (SHOW_FINDER(hwin->role)) {
	vwin_add_finder(hwin);
    }

    gtk_widget_show_all(hbox);
}

void show_gui_help (int helpcode)
{
    int idx = gui_ci_to_index(helpcode);
    int pos, role = GUI_HELP;

    /* try for GUI help first */
    pos = help_pos_from_index(idx, role);

    if (pos <= 0 && translated_cmdref) {
	/* English GUI help? */
	role = GUI_HELP_EN;
	pos = help_pos_from_index(idx, role);
    }

    if (pos <= 0) {
	/* command help? */
	role = CMD_HELP;
	pos = help_pos_from_index(idx, role);
    }

    if (pos <= 0 && translated_cmdref) {
	/* English command help? */
	role = CMD_HELP_EN;
	pos = help_pos_from_index(idx, role);
    }

    if (pos > 0) {
	real_do_help(idx, pos, role);
    } else {
	warnbox(_("Sorry, no help is available"));
    }
}

void context_help (GtkWidget *widget, gpointer data)
{
    int helpcode = GPOINTER_TO_INT(data);

    show_gui_help(helpcode);
}

GtkWidget *context_help_button (GtkWidget *hbox, int helpcode)
{
    GtkWidget *button;

    button = gtk_button_new_from_stock(GTK_STOCK_HELP);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       button, TRUE);

    if (helpcode >= 0) {
	/* passing helpcode < 0 is a way of reserving
	   the callback for something special
	*/
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(context_help),
			 GINT_TO_POINTER(helpcode));
    }

    return button;
}

static gboolean nullify_hwin (GtkWidget *w, windata_t **phwin)
{
    *phwin = NULL;
    return FALSE;
}

/* sync the tree index view with the currently selected topic, if it's
   not already in sync */

static void helpwin_set_topic_index (windata_t *hwin, int idx)
{
    GtkWidget *w =
	g_object_get_data(G_OBJECT(hwin->text), "tview");
    GtkTreeView *view = GTK_TREE_VIEW(w);

    hwin->active_var = idx;

    if (view != NULL) {
	GtkTreeModel *model = gtk_tree_view_get_model(view);
	GtkTreeIter iter, parent;
	gboolean ok;

	if (idx == 0) {
	    ok = gtk_tree_model_get_iter_first(model, &iter);
	} else {
	    ok = help_iter_from_index(idx, hwin->role, &iter,
				      &parent);
	}

	if (ok) {
	    GtkTreeSelection *sel;

	    sel = gtk_tree_view_get_selection(view);
	    if (!gtk_tree_selection_iter_is_selected(sel, &iter)) {
		GtkTreePath *path;

		/* gtk_tree_view_collapse_all(view); should we? */
		gtk_tree_selection_select_iter(sel, &iter);
		path = gtk_tree_model_get_path(model, &iter);
		gtk_tree_view_expand_to_path(view, path);
		gtk_tree_view_set_cursor(view, path, NULL, FALSE);
		gtk_tree_path_free(path);
	    }
	}
    }
}

static windata_t *real_do_help (int idx, int pos, int role)
{
    static windata_t *gui_hwin;
    static windata_t *cmd_hwin;
    static windata_t *func_hwin;
    static windata_t *en_gui_hwin;
    static windata_t *en_cmd_hwin;
    static windata_t *en_func_hwin;
    windata_t *hwin = NULL;
    const char *fname = NULL;

    if (pos < 0) {
	dummy_call();
	return NULL;
    }

#if HDEBUG
    fprintf(stderr, "real_do_help: idx=%d, pos=%d, role=%d\n",
	    idx, pos, role);
    fprintf(stderr, "gui_hwin = %p\n", (void *) gui_hwin);
    fprintf(stderr, "cmd_hwin = %p\n", (void *) cmd_hwin);
    fprintf(stderr, "func_hwin = %p\n", (void *) func_hwin);
    fprintf(stderr, "en_gui_hwin = %p\n", (void *) en_gui_hwin);
    fprintf(stderr, "en_cmd_hwin = %p\n", (void *) en_cmd_hwin);
    fprintf(stderr, "en_func_hwin = %p\n", (void *) en_func_hwin);
#endif

    if (role == CMD_HELP) {
	hwin = cmd_hwin;
	fname = helpfile_path(GRETL_CMDREF, 0, 0);
    } else if (role == GUI_HELP) {
	hwin = gui_hwin;
	fname = helpfile_path(GRETL_GUI_HELP, 0, 0);
    } else if (role == FUNC_HELP) {
	hwin = func_hwin;
	fname = helpfile_path(GRETL_FUNCREF, 0, 0);
    } else if (role == CMD_HELP_EN) {
	hwin = en_cmd_hwin;
	fname = helpfile_path(GRETL_CMDREF, 0, 1);
    } else if (role == GUI_HELP_EN) {
	hwin = en_gui_hwin;
	fname = helpfile_path(GRETL_GUI_HELP, 0, 1);
    } else if (role == FUNC_HELP_EN) {
	hwin = en_func_hwin;
	fname = helpfile_path(GRETL_FUNCREF, 0, 1);
    }

    if (hwin != NULL) {
	gtk_window_present(GTK_WINDOW(hwin->main));
    } else {
	hwin = view_help_file(fname, role);
	if (hwin != NULL) {
	    windata_t **phwin = NULL;

	    if (role == CMD_HELP) {
		cmd_hwin = hwin;
		phwin = &cmd_hwin;
	    } else if (role == GUI_HELP) {
		gui_hwin = hwin;
		phwin = &gui_hwin;
	    } else if (role == FUNC_HELP) {
		func_hwin = hwin;
		phwin = &func_hwin;
	    } else if (role == CMD_HELP_EN) {
		en_cmd_hwin = hwin;
		phwin = &en_cmd_hwin;
	    } else if (role == GUI_HELP_EN) {
		en_gui_hwin = hwin;
		phwin = &en_gui_hwin;
	    } else if (role == FUNC_HELP_EN) {
		en_func_hwin = hwin;
		phwin = &en_func_hwin;
	    }

	    g_signal_connect(G_OBJECT(hwin->main), "destroy",
			     G_CALLBACK(nullify_hwin), phwin);
	}
    }

#if HDEBUG
    fprintf(stderr, "real_do_help: doing set_help_topic_buffer:\n"
	    " hwin=%p, pos=%d, role=%d\n", (void *) hwin, pos, role);
#endif

    if (hwin != NULL) {
	int ret = set_help_topic_buffer(hwin, pos);

	if (ret >= 0) {
	    helpwin_set_topic_index(hwin, idx);
	}
    }

    return hwin;
}

void display_text_help (GtkAction *action)
{
    if (action != NULL) {
	const char *aname = gtk_action_get_name(action);

	if (!strcmp(aname, "TextCmdRef")) {
	    real_do_help(0, 0, CMD_HELP);
	} else if (!strcmp(aname, "FuncRef")) {
	    real_do_help(0, 0, FUNC_HELP);
	} else if (!strcmp(aname, "PkgHelp")) {
	    show_gui_help(PKGHELP);
	}
    } else {
	/* default: commands help */
	real_do_help(0, 0, CMD_HELP);
    }
}

/* called from textbuf.c */

void command_help_callback (int idx, int en)
{
    int role = (en)? CMD_HELP_EN : CMD_HELP;
    int pos = 0;

    if (idx > NC) {
	/* a GUI-special help item */
	show_gui_help(idx);
	return;
    }

    if (idx > 0) {
	pos = help_pos_from_index(idx, role);
	if (pos < 0 && !en && translated_cmdref) {
	    /* no translated entry: fall back on English */
	    role = CMD_HELP_EN;
	    pos = help_pos_from_index(idx, role);
	}
    }

    real_do_help(idx, pos, role);
}

/* called from textbuf.c */

void function_help_callback (int idx, int en)
{
    int role = (en)? FUNC_HELP_EN : FUNC_HELP;
    int pos = help_pos_from_index(idx, role);

    real_do_help(idx, pos, role);
}

/* below: must return > 0 to do anything useful */

static int help_pos_from_string (const char *s, int *idx, int *role)
{
    char word[16];
    int pos = 0;

    *word = '\0';
    strncat(word, s, 15);

    if (*role != FUNC_HELP) {
	*idx = gretl_command_number(word);
	pos = help_pos_from_index(*idx, *role);
	if (pos <= 0 && translated_cmdref) {
	    pos = help_pos_from_index(*idx, CMD_HELP_EN);
	    if (pos > 0) {
		*role = CMD_HELP_EN;
	    }
	}
    }

    if (pos <= 0) {
	/* try for function instead of command? */
	pos = function_help_pos_from_word(word, FUNC_HELP);
	if (pos > 0) {
	    *role = FUNC_HELP;
	    *idx = function_help_index_from_word(word, *role);
	} else if (translated_fnref) {
	    pos = function_help_pos_from_word(word, FUNC_HELP_EN);
	    if (pos > 0) {
		*role = FUNC_HELP_EN;
		*idx = function_help_index_from_word(word, *role);
	    }
	}
    }

    return pos;
}

/* Having found a 'word' at the cursor, see whether it's immediately
   preceded by '$': if so, the word should be extended backward to
   include that character.
*/

static int is_dollar_word (GtkTextBuffer *buf,
			   GtkTextIter *w_start,
			   gchar **textp)
{
    GtkTextIter d_start = *w_start;
    int ret = 0;

    if (gtk_text_iter_backward_char(&d_start)) {
	gchar *dtest = gtk_text_buffer_get_text(buf, &d_start,
						w_start, FALSE);

	if (*dtest == '$') {
	    gchar *s = g_strdup_printf("$%s", *textp);

	    g_free(*textp);
	    *textp = s;
	    ret = 1;
	}
	g_free(dtest);
    }

    return ret;
}

/* In case a given identifier has both a command and a function form,
   try to determine if we're looking at the function form, using the
   fact that the identifier can be followed by left parenthesis only
   in the function case.
*/

static int probably_function (GtkTextBuffer *buf, GtkTextIter *w_end)
{
    gunichar gu = gtk_text_iter_get_char(w_end);
    gchar *chk = g_ucs4_to_utf8(&gu, 1, NULL, NULL, NULL);
    int ret = 0;

    if (chk != NULL) {
        ret = *chk == '(';
	g_free(chk);
    }

    return ret;
}

gint interactive_script_help (GtkWidget *widget, GdkEventButton *b,
			      windata_t *vwin)
{
    if (!window_help_is_active(vwin)) {
	/* command help not activated */
	return FALSE;
    } else {
	gchar *text = NULL;
	int role = CMD_HELP;
	GtkTextBuffer *buf;
	GtkTextIter iter;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
	gtk_text_buffer_get_iter_at_mark(buf, &iter,
					 gtk_text_buffer_get_insert(buf));

	if (gtk_text_iter_inside_word(&iter)) {
	    GtkTextIter w_start = iter;
	    GtkTextIter w_end = iter;

	    if (!gtk_text_iter_starts_word(&iter)) {
		gtk_text_iter_backward_word_start(&w_start);
	    }
	    if (!gtk_text_iter_ends_word(&iter)) {
		gtk_text_iter_forward_word_end(&w_end);
	    }
	    text = gtk_text_buffer_get_text(buf, &w_start, &w_end, FALSE);

	    if (text != NULL) {
		if (is_dollar_word(buf, &w_start, &text)) {
		    /* got $<word> */
		    role = FUNC_HELP;
		} else if (probably_function(buf, &w_end)) {
		    /* got <word>( */
		    role = FUNC_HELP;
		}
	    }
	}

	if (text != NULL && *text != '\0') {
	    int idx, pos;

	    pos = help_pos_from_string(text, &idx, &role);
	    if (pos <= 0) {
		warnbox(_("Sorry, help not found"));
	    } else {
		real_do_help(idx, pos, role);
	    }
	}

	/* clean up */
	unset_window_help_active(vwin);
	text_set_cursor(vwin->text, 0);
	g_free(text);
    }

    return FALSE;
}

/* First response to "help <param>" in GUI console, when given with no
   option: if we got a command word or function name, pop open a
   nicely formatted help window.  If this function returns non-zero
   we'll fall back on the command-line help function.
*/

int gui_console_help (const char *param)
{
    int idx = 0, role = CMD_HELP;
    int pos, err = 0;

    pos = help_pos_from_string(param, &idx, &role);

    if (pos <= 0) {
	err = 1;
    } else {
	real_do_help(idx, pos, role);
    }

    return err;
}

void text_find (gpointer unused, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin->finder != NULL) {
	GtkWidget *p = gtk_widget_get_parent(vwin->finder);

	if (!gtk_widget_get_visible(p)) {
	    gtk_widget_show_all(p);
	}
	gtk_widget_grab_focus(vwin->finder);
	gtk_editable_select_region(GTK_EDITABLE(vwin->finder),
				   0, -1);
    } else if (vwin->flags & VWIN_USE_FOOTER) {
	vwin_add_footer_finder(vwin);
    } else {
	find_string_dialog(find_in_text, vwin);
    }
}

void text_find_again (gpointer unused, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin->finder != NULL) {
	if (gtk_widget_get_visible(vwin->finder)) {
	    g_signal_emit_by_name(G_OBJECT(vwin->finder), "activate", NULL);
	}
    } else if (find_dialog != NULL) {
	if (vwin == g_object_get_data(G_OBJECT(find_dialog), "windat")) {
	    find_in_text(NULL, find_dialog);
	}
    }
}

#ifndef GRETL_EDIT

void listbox_find (gpointer unused, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin->finder != NULL) {
	gtk_widget_grab_focus(vwin->finder);
	gtk_editable_select_region(GTK_EDITABLE(vwin->finder),
				   0, -1);
    } else {
	find_string_dialog(find_in_listbox, vwin);
    }
}

void listbox_find_again (gpointer unused, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin->finder != NULL) {
	g_signal_emit_by_name(G_OBJECT(vwin->finder), "activate", NULL);
    } else if (find_dialog != NULL) {
	if (vwin == g_object_get_data(G_OBJECT(find_dialog), "windat")) {
	    find_in_listbox(NULL, find_dialog);
	}
    }
}

static int string_match_pos (const char *haystack, const char *needle,
			     gboolean sensitive, int start)
{
    int hlen = strlen(haystack);
    int nlen = strlen(needle);
    int pos, found;

    for (pos = start; pos < hlen; pos++) {
	if (sensitive) {
	    found = !strncmp(&haystack[pos], needle, nlen);
	} else {
	    found = !g_ascii_strncasecmp(&haystack[pos], needle, nlen);
	}
        if (found) {
             return pos;
	}
    }

    return -1;
}

#endif /* not GRETL_EDIT */

static gint close_find_dialog (GtkWidget *widget, gpointer data)
{
    find_dialog = NULL;
    return FALSE;
}

/* search for string @s in text buffer associated with @view */

static gboolean real_find_in_text (GtkTextView *view, const gchar *s,
				   gboolean sensitive,
				   gboolean from_cursor,
				   gboolean search_all)
{
    GtkTextBuffer *buf;
    GtkTextIter iter, start, end;
    GtkTextMark *vis;
    int found = 0;
    int wrapped = 0;
    int n = strlen(s);
    gchar *got;

    buf = gtk_text_view_get_buffer(view);

 text_search_wrap:

    while (gtk_events_pending()) {
	gtk_main_iteration();
    }

    if (from_cursor) {
	GtkTextIter sel_bound;

	gtk_text_buffer_get_iter_at_mark(buf, &iter,
					 gtk_text_buffer_get_insert(buf));
	gtk_text_buffer_get_iter_at_mark(buf, &sel_bound,
					 gtk_text_buffer_get_selection_bound(buf));
	gtk_text_iter_order(&sel_bound, &iter);
    } else {
	gtk_text_buffer_get_iter_at_offset(buf, &iter, 0);
    }

    start = end = iter;

    if (!gtk_text_iter_forward_chars(&end, n)) {
	/* we're already at end of the buffer */
	if (from_cursor && !wrapped && !search_all) {
	    from_cursor = FALSE;
	    wrapped = 1;
	    goto text_search_wrap;
	} else {
	    return 0;
	}
    }

    while (!found) {
	got = gtk_text_buffer_get_text(buf, &start, &end, FALSE);
	if (sensitive) {
	    found = !strcmp(got, s);
	} else {
	    found = !g_ascii_strcasecmp(got, s);
	}
	g_free(got);
	if (found || !gtk_text_iter_forward_char(&start) ||
	    !gtk_text_iter_forward_char(&end)) {
	    break;
	}
    }

    if (found) {
	gtk_text_buffer_place_cursor(buf, &start);
	gtk_text_buffer_move_mark_by_name(buf, "selection_bound", &end);
	vis = gtk_text_buffer_create_mark(buf, "vis", &end, FALSE);
	gtk_text_view_scroll_to_mark(view, vis, 0.05, FALSE, 0, 0);
    } else if (from_cursor && !wrapped && !search_all) {
	/* try wrapping */
	from_cursor = FALSE;
	wrapped = 1;
	goto text_search_wrap;
    }

    return found;
}

static void find_in_text (GtkWidget *button, GtkWidget *dialog)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(dialog), "windat");
    gboolean found, sensitive;

    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);
    if (needle == NULL || *needle == '\0') {
	return;
    }

    sensitive = !all_lower_case(needle);

    found = real_find_in_text(GTK_TEXT_VIEW(vwin->text), needle,
			      sensitive, TRUE, FALSE);

    if (!found) {
	notify_string_not_found(find_entry);
    }
}

#ifndef GRETL_EDIT

static gboolean real_find_in_listbox (windata_t *vwin,
				      const gchar *s,
				      gboolean sensitive,
				      gboolean vnames)
{
    int search_cols[4] = {0, 0, -1, -1};
    int minvar, wrapped = 0;
    gchar *haystack;
    char pstr[16];
    GtkTreeModel *model = NULL;
    GtkTreeIter iter;
    gboolean got_iter;
    int i, pos = -1;

    /* first check that there's something to search */
    if (vwin->listbox != NULL) {
	model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
    }
    if (model == NULL) {
	return FALSE;
    }

    /* if searching in the main gretl window, start on line 1 */
    minvar = (vwin == mdata)? 1 : 0;

    /* first try to get the current line plus one as starting point */
    sprintf(pstr, "%d", vwin->active_var);
    got_iter = gtk_tree_model_get_iter_from_string(model, &iter, pstr);
    if (got_iter) {
	got_iter = gtk_tree_model_iter_next(model, &iter);
    }

    if (!got_iter) {
	/* fallback: start from the top */
	got_iter = gtk_tree_model_get_iter_first(model, &iter);
    }

    if (!got_iter) {
	/* failed totally, get out */
	return FALSE;
    }

    if (vnames) {
	/* case-sensitive search for series names */
	search_cols[0] = 1;  /* series name */
	search_cols[1] = -1; /* invalid */
    } else if (vwin == mdata) {
	search_cols[0] = 1;  /* series name */
	search_cols[1] = 2;  /* description */
    } else if (vwin->role == FUNC_FILES) {
	search_cols[0] = 0;  /* package name */
	search_cols[1] = 4;  /* description */
	search_cols[2] = 3;  /* author */
    } else if (vwin->role == REMOTE_FUNC_FILES) {
	search_cols[0] = 0;  /* package name */
	search_cols[1] = 4;  /* description */
	search_cols[2] = 3;  /* author */
    } else if (vwin->role == DBNOMICS_DB) {
	search_cols[0] = 1;  /* content */
	search_cols[1] = 0;  /* code */
    } else {
	/* databases, datafiles */
	search_cols[0] = 1; /* description */
	search_cols[1] = 0; /* filename */
    }

 search_wrap:

    while (pos < 0) {
	for (i=0; pos < 0 && search_cols[i] >= 0; i++) {
	    /* iterate across searchable columns */
	    gtk_tree_model_get(model, &iter, search_cols[i], &haystack, -1);
	    if (haystack != NULL) {
		if (*haystack != '\0') {
		    pos = string_match_pos(haystack, needle, sensitive, 0);
		}
		g_free(haystack);
		haystack = NULL;
	    }
	}
	if (pos >= 0 || !gtk_tree_model_iter_next(model, &iter)) {
	    /* either found, or nowhere left to search */
	    break;
	}
    }

    if (pos < 0 && vwin->active_var > minvar && !wrapped) {
	/* try wrapping to start */
	gtk_tree_model_get_iter_first(model, &iter);
	if (minvar > 0 && !gtk_tree_model_iter_next(model, &iter)) {
	    ; /* do nothing: there's only one line in the box */
	} else {
	    wrapped = 1;
	    goto search_wrap;
	}
    }

    if (pos >= 0) {
	GtkTreePath *path = gtk_tree_model_get_path(model, &iter);

	gtk_tree_view_set_cursor(GTK_TREE_VIEW(vwin->listbox),
				 path, NULL, FALSE);
	gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(vwin->listbox),
				     path, NULL, TRUE, 0.5, 0);
	vwin->active_var = tree_path_get_row_number(path);
	gtk_tree_path_free(path);
    }

    return (pos >= 0);
}

/* Given @targ, the name of a function package, try to find
   it in @vwin's listbox, and if found focus that row.
   Return TRUE if found, FALSE otherwise.
*/

gboolean find_package_in_viewer (windata_t *vwin,
				 const gchar *targ)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *haystack;
    int pos = -1;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));

    if (model == NULL || !gtk_tree_model_get_iter_first(model, &iter)) {
	return FALSE;
    }

    while (pos < 0) {
	gtk_tree_model_get(model, &iter, 0, &haystack, -1);
	if (haystack != NULL) {
	    if (*haystack != '\0') {
		pos = string_match_pos(haystack, targ, TRUE, 0);
	    }
	    g_free(haystack);
	    haystack = NULL;
	    if (pos >= 0 || !gtk_tree_model_iter_next(model, &iter)) {
		break;
	    }
	}
    }

    if (pos >= 0) {
	GtkTreePath *path = gtk_tree_model_get_path(model, &iter);

	gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(vwin->listbox),
				     path, NULL, FALSE, 0, 0);
	gtk_tree_view_set_cursor(GTK_TREE_VIEW(vwin->listbox),
				 path, NULL, FALSE);
	vwin->active_var = tree_path_get_row_number(path);
	gtk_tree_path_free(path);
    }

    return (pos >= 0);
}

/* used for windows that do not have a built-in search entry,
   but which call the function find_string_dialog() */

static void find_in_listbox (GtkWidget *w, GtkWidget *dialog)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(dialog), "windat");
    gpointer vp;
    gboolean sensitive;
    gboolean vnames = FALSE;
    gboolean found;

    if (needle != NULL) {
	g_free(needle);
	needle = NULL;
    }

    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);
    if (needle == NULL || *needle == '\0') {
	return;
    }

    sensitive = !all_lower_case(needle);

    /* are we confining the search to variable names? */
    vp = g_object_get_data(G_OBJECT(dialog), "vnames_only");
    if (vp != NULL) {
	vnames = GPOINTER_TO_INT(vp);
	if (vnames) {
	    /* varname search is advertised as case-sensitive */
	    sensitive = TRUE;
	}
    }

    found = real_find_in_listbox(vwin, needle, sensitive, vnames);

    if (!found) {
	notify_string_not_found(find_entry);
    }
}

#endif /* not GRETL_EDIT */

static void cancel_find (GtkWidget *button, GtkWidget *dialog)
{
    if (find_dialog != NULL) {
	gtk_widget_destroy(dialog);
	find_dialog = NULL;
    }
}

static void parent_find (GtkWidget *finder, windata_t *caller)
{
    GtkWidget *w = vwin_toplevel(caller);

    if (w != NULL) {
	gtk_window_set_transient_for(GTK_WINDOW(finder), GTK_WINDOW(w));
	gtk_window_set_destroy_with_parent(GTK_WINDOW(finder), TRUE);
    }
}

static void toggle_vname_search (GtkToggleButton *tb, GtkWidget *w)
{
    if (gtk_toggle_button_get_active(tb)) {
	g_object_set_data(G_OBJECT(w), "vnames_only",
			  GINT_TO_POINTER(1));
    } else {
	g_object_set_data(G_OBJECT(w), "vnames_only",
			  GINT_TO_POINTER(0));
    }
}

static gint maybe_find_again (GtkWidget *w, GdkEventKey *event,
			      GtkWidget *button)
{
    if ((event->state & GDK_CONTROL_MASK) &&
	(event->keyval == GDK_g || event->keyval == GDK_G)) {
	g_signal_emit_by_name(G_OBJECT(button), "clicked", NULL);
	return TRUE;
    } else {
	return FALSE;
    }
}

static void find_string_dialog (void (*findfunc)(), windata_t *vwin)
{
    GtkWidget *parent;
    GtkWidget *label;
    GtkWidget *button;
    GtkWidget *vbox;
    GtkWidget *hbox;

    if (find_dialog != NULL) {
	g_object_set_data(G_OBJECT(find_dialog), "windat", vwin);
	parent_find(find_dialog, vwin);
	gtk_window_present(GTK_WINDOW(find_dialog));
	return;
    }

    parent = vwin->topmain != NULL ? vwin->topmain : vwin->main;
    find_dialog = gretl_dialog_new(_("gretl: find"), parent, 0);
    g_object_set_data(G_OBJECT(find_dialog), "windat", vwin);

    g_signal_connect(G_OBJECT(find_dialog), "destroy",
		     G_CALLBACK(close_find_dialog),
		     find_dialog);

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_(" Find what:"));
    gtk_widget_show(label);
    find_entry = gtk_entry_new();

    if (needle != NULL) {
	gtk_entry_set_text(GTK_ENTRY(find_entry), needle);
	gtk_editable_select_region(GTK_EDITABLE(find_entry), 0, -1);
    }

    g_signal_connect(G_OBJECT(find_entry), "activate",
		     G_CALLBACK(findfunc), find_dialog);
    gtk_widget_show(find_entry);
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), find_entry, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(find_dialog));

    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    if (vwin == mdata) {
	hbox = gtk_hbox_new(FALSE, 5);
	button = gtk_check_button_new_with_label(_("Variable names only (case sensitive)"));
	g_signal_connect(G_OBJECT(button), "toggled",
			 G_CALLBACK(toggle_vname_search), find_dialog);
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
	gtk_widget_show_all(hbox);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(find_dialog));

    /* Close button */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(cancel_find), find_dialog);
    gtk_widget_show(button);

    /* Find button */
    button = gtk_button_new_from_stock(GTK_STOCK_FIND);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(findfunc), find_dialog);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    g_signal_connect(G_OBJECT(find_entry), "key-press-event",
		     G_CALLBACK(maybe_find_again), button);

    gtk_widget_grab_focus(find_entry);
    gtk_widget_show(find_dialog);
}

enum {
    GRETL_GUIDE = 1,
    GRETL_REF,
    GNUPLOT_REF,
    X12A_REF,
    GRETL_KEYS,
    HANSL_PRIMER,
    PKGBOOK,
    GRETL_MPI,
    GRETL_SVM,
    GRETL_DBN,
    GRETL_GEO,
    GRETL_LP
};

static int get_writable_doc_path (char *path, const char *fname)
{
    static int sysdoc_writable = -1;
    static int userdoc_writable = -1;
    const char *gretldir = gretl_home();
    const char *dotdir = gretl_dotdir();
    FILE *fp;
    int err = 0;

#ifdef G_OS_WIN32
    sysdoc_writable = 0;
#endif

    if (sysdoc_writable == 1) {
	sprintf(path, "%sdoc%c%s", gretldir, SLASH, fname);
	return 0;
    } else if (userdoc_writable == 1) {
	sprintf(path, "%sdoc%c%s", dotdir, SLASH, fname);
	return 0;
    }

    if (sysdoc_writable < 0) {
	sysdoc_writable = 0;
	sprintf(path, "%sdoc", gretldir);
	if (gretl_mkdir(path) == 0) {
	    strcat(path, SLASHSTR);
	    strcat(path, fname);
	    fp = gretl_fopen(path, "w");
	    if (fp != NULL) {
		sysdoc_writable = 1;
		fclose(fp);
		gretl_remove(path);
	    }
	}
    }

    if (!sysdoc_writable && userdoc_writable < 0) {
	/* can't write to 'sys' dir, user dir not tested yet */
	userdoc_writable = 0;
	sprintf(path, "%sdoc", dotdir);
	if (gretl_mkdir(path) == 0) {
	    sprintf(path, "%sdoc%c%s", dotdir, SLASH, fname);
	    fp = gretl_fopen(path, "w");
	    if (fp != NULL) {
		userdoc_writable = 1;
		fclose(fp);
		gretl_remove(path);
	    }
	}
    }

    if (!sysdoc_writable && !userdoc_writable) {
	err = 1;
    }

    return err;
}

static int get_x12a_doc_path (char *path, const char *fname)
{
    const char *x12a = gretl_x12_arima();
    int ret = 0;

    *path = '\0';

    if (x12a != NULL && *x12a != '\0') {
	char *p;

	strcpy(path, x12a);
	p = strrslash(path);
	if (p != NULL) {
	    sprintf(p + 1, "docs%c%s", SLASH, fname);
	    ret = 1;
	} else {
	    *path = '\0';
	}

#if !defined(G_OS_WIN32) && !defined(OS_OSX)
	if (!ret) {
	    /* using gretl x12a package? */
	    if (gretl_x12_is_x13()) {
		sprintf(path, "/opt/x13as/docs/%s", fname);
	    } else {
		sprintf(path, "/opt/x12arima/docs/%s", fname);
	    }
	    ret = 1;
	}
#endif
    }

#if !defined(G_OS_WIN32)
    if (ret && gretl_test_fopen(path, "r") != 0) {
	/* try lower-casing filename (recent x13as!) */
	char *tmp = gretl_strdup(fname);

	sprintf(path, "/opt/x13as/docs/%s", gretl_lower(tmp));
	free(tmp);
    }
#endif

    return ret;
}

/* Get a language-specific query string, for asking
   the user whether a translation is preferred.
*/

static const char *tr_query (const char *lang)
{
    if (!strcmp(lang, "es")) {
	return "¿Mostrar traducción al español?";
    } else if (!strcmp(lang, "gl")) {
	return "Mostrar tradución ao galego?";
    } else if (!strcmp(lang, "it")) {
	return "Mostra traduzione in italiano?";
    } else if (!strcmp(lang, "pt")) {
	return "Mostrar portugues tradução?";
    } else if (!strcmp(lang, "ru")) {
	return "Показать перевод на русский язык?";
    } else {
	return NULL;
    }
}

/* For @code giving the ID number of a doc resource and
   @lang identifying a language, return the filename of a
   language-specific version of the resource, or NULL if
   none is available.
*/

static const char *have_translation (int code, const char *lang)
{
    gchar *ret = NULL;

    if (code == HANSL_PRIMER) {
	if (!strcmp(lang, "ru")) {
	    ret = "hansl-primer-ru.pdf";
	}
    } else if (code == GRETL_REF) {
	if (!strcmp(lang, "es")) {
	    ret = "gretl-ref-es.pdf";
	} else if (!strcmp(lang, "gl")) {
	    ret = "gretl-ref-gl.pdf";
	} else if (!strcmp(lang, "it")) {
	    ret = "gretl-ref-it.pdf";
	} else if (!strcmp(lang, "pt")) {
	    ret = "gretl-ref-pt.pdf";
	}
    }

    return ret;
}

/* Determine if we should show a translation of the
   doc resource identified by @code. If so, return
   the required filename; if not, return NULL.
*/

static const char *show_translation (int code)
{
    const char *fname = NULL;
    char lang[3] = {0};

#ifdef WIN32
    gchar *loc = g_win32_getlocale();

    strncat(lang, loc, 2);
    if (loc != NULL) {
	fname = have_translation(code, lang);
	g_free(loc);
    }
#elif defined(ENABLE_NLS)
    char *loc = setlocale(LC_MESSAGES, NULL);

    if (loc != NULL) {
	strncat(lang, loc, 2);
	fname = have_translation(code, lang);
    }
#endif

    if (fname != NULL) {
	/* We have a translation, but does the user want it? */
	const char *msg = tr_query(lang);
	int resp = yes_no_dialog(NULL, msg, NULL);

	if (resp != GRETL_YES) {
	    fname = NULL;
	}
    }

    return fname;
}

/* @pref is the documentation preference registered in settings.c:
   0 = English, US letter
   1 = English, A4
   [2 = Translation, if available]
*/

static int find_or_download_pdf (int code, int pref, char *fullpath)
{
    const char *guide_files[] = {
	"gretl-guide.pdf",
	"gretl-guide-a4.pdf"
    };
    const char *ref_files[] = {
	"gretl-ref.pdf",
	"gretl-ref-a4.pdf",
    };
    const char *kbd_files[] = {
	"gretl-keys.pdf",
	"gretl-keys-a4.pdf"
    };
    const char *primer_files[] = {
	"hansl-primer.pdf",
	"hansl-primer-a4.pdf",
    };
    const char *pkgbook_files[] = {
	"pkgbook.pdf",
	"pkgbook-a4.pdf"
    };
    const char *gretlMPI_files[] = {
	"gretl-mpi.pdf",
	"gretl-mpi-a4.pdf"
    };
    const char *gretlSVM_files[] = {
	"gretl-svm.pdf",
	"gretl-svm-a4.pdf"
    };
    const char *gretlLP_files[] = {
	"gretl-lpsolve.pdf",
	"gretl-lpsolve-a4.pdf"
    };
    const char *fname = NULL;
    int gotit = 0;
    int err = 0;

    if (pref < 0 || pref > 2) {
	/* out of bounds */
	pref = 0;
    }

    if (pref > 0) {
	/* Try offering a translation where available: currently only
	   for the Gretl Reference and Hansl primer (Russian).
	*/
	pref = 1;
	if (code == HANSL_PRIMER || code == GRETL_REF) {
	    fname = show_translation(code);
	}
    }

#if 0
    fprintf(stderr, "HERE code=%d, pref=%d, fname %s\n",
	    code, pref, fname != NULL ? fname : "TBD");
#endif

    if (fname != NULL) {
	/* we got a specific translation */
	goto next_step;
    }

    if (code == GRETL_GUIDE) {
	fname = guide_files[pref];
    } else if (code == GRETL_REF) {
	fname = ref_files[pref];
    } else if (code == GRETL_KEYS) {
	fname = kbd_files[pref];
    } else if (code == HANSL_PRIMER) {
	fname = primer_files[pref];
    } else if (code == PKGBOOK) {
	fname = pkgbook_files[pref];
    } else if (code == GRETL_MPI) {
	fname = gretlMPI_files[pref];
    } else if (code == GRETL_SVM) {
	fname = gretlSVM_files[pref];
    } else if (code == GRETL_LP) {
	fname = gretlLP_files[pref];
    } else if (code == GNUPLOT_REF) {
	fname = "gnuplot.pdf";
    } else if (code == X12A_REF) {
	fname = gretl_x12_is_x13() ? "docX13AS.pdf" : "x12adocV03.pdf";
    } else if (code == GRETL_DBN) {
	fname = "dbnomics.pdf";
	sprintf(fullpath, "%sfunctions%cdbnomics%c%s",
		gretl_home(), SLASH, SLASH, fname);
    } else if (code == GRETL_GEO) {
	fname = "geoplot.pdf";
	sprintf(fullpath, "%sfunctions%cgeoplot%c%s",
		gretl_home(), SLASH, SLASH, fname);
    } else {
	return E_DATA;
    }

 next_step:

    fprintf(stderr, "pdf help: looking for '%s'\n", fname);

    if (code != GRETL_DBN && code != GRETL_GEO) {
	/* is the file available in public dir? */
	sprintf(fullpath, "%sdoc%c%s", gretl_home(), SLASH, fname);
    }

    if (*fullpath != '\0' && gretl_test_fopen(fullpath, "r") == 0) {
	gotit = 1;
    } else if (fname != NULL && *fname != '\0' &&
	       gretl_test_fopen(fname, "r") == 0) {
	gotit = 1;
    }

    if (!gotit && code == X12A_REF) {
	*fullpath = '\0';
	get_x12a_doc_path(fullpath, fname);
	if (*fullpath != '\0') {
	    err = gretl_test_fopen(fullpath, "r");
	    if (!err) {
		gotit = 1;
	    }
	}
    }

    if (!gotit) {
	/* try in the user's dotdir? */
	if (code == GRETL_DBN) {
	    sprintf(fullpath, "%sfunctions%cdbnomics%cdbnomics.pdf",
		    gretl_dotdir(), SLASH, SLASH);
	} else if (code == GRETL_GEO) {
	    sprintf(fullpath, "%sfunctions%cgeojson%cgeoplot.pdf",
		    gretl_dotdir(), SLASH, SLASH);
	} else {
	    sprintf(fullpath, "%sdoc%c%s", gretl_dotdir(), SLASH, fname);
	}
	err = gretl_test_fopen(fullpath, "r");
	if (!err) {
	    gotit = 1;
	}
    }

#ifndef GRETL_EDIT
    if (!gotit && code == GRETL_DBN) {
	/* try installing the dbnomics package */
	char *dlpath = NULL;

	err = download_addon("dbnomics", &dlpath);
	if (!err) {
	    /* .gfn -> .pdf */
	    switch_ext(fullpath, dlpath, "pdf");
	    free(dlpath);
	}
    }
#endif
    if (!gotit && code != GRETL_DBN) {
	/* try downloading the manual file */
	err = get_writable_doc_path(fullpath, fname);
	if (!err) {
	    err = retrieve_manfile(fname, fullpath);
	}
    }

    if (err) {
	const char *buf = gretl_errmsg_get();

	if (*buf) {
	    errbox(buf);
	} else {
	    errbox(_("Failed to download file"));
	}
    }

    return err;
}

int get_pdf_path (const char *name, char *fullpath)
{
    int code = 0;

    if (!strcmp(name, "gretl-lpsolve.pdf")) {
	code = GRETL_LP;
    } else if (!strcmp(name, "gretl-svm.pdf")) {
	code = GRETL_SVM;
    } else if (!strcmp(name, "gretl-mpi.pdf")) {
	code = GRETL_MPI;
    }

    if (code > 0) {
	return find_or_download_pdf(code, 0, fullpath);
    } else {
	return 1;
    }
}

void gretl_show_pdf (const char *fname, const char *option)
{
#if defined(G_OS_WIN32)
    if (option != NULL) {
	win32_open_pdf(fname, option);
    } else {
	win32_open_file(fname);
    }
#elif defined(OS_OSX)
    if (option != NULL) {
	osx_open_pdf(fname, option);
    } else {
	osx_open_file(fname);
    }
#else
    gretl_fork("viewpdf", fname, option);
#endif
}

void display_pdf_help (GtkAction *action)
{
    char fname[FILENAME_MAX];
    int err, code = GRETL_GUIDE;

    if (action != NULL) {
	const char *aname = gtk_action_get_name(action);

	if (!strcmp(aname, "PDFCmdRef")) {
	    code = GRETL_REF;
	} else if (!strcmp(aname, "KbdRef")) {
	    code = GRETL_KEYS;
	} else if (!strcmp(aname, "Primer")) {
	    code = HANSL_PRIMER;
	} else if (!strcmp(aname, "Pkgbook")) {
	    code = PKGBOOK;
	} else if (!strcmp(aname, "GeoplotDoc")) {
	    code = GRETL_GEO;
	} else if (!strcmp(aname, "gretlMPI")) {
	    code = GRETL_MPI;
	} else if (!strcmp(aname, "gretlSVM")) {
	    code = GRETL_SVM;
	} else if (!strcmp(aname, "gretlDBN")) {
	    code = GRETL_DBN;
	} else if (!strcmp(aname, "gretlLpsolve")) {
	    code = GRETL_LP;
	}
    }

    err = find_or_download_pdf(code, get_manpref(), fname);

    if (!err) {
	gretl_show_pdf(fname, NULL);
    }
}

void display_guide_chapter (const char *dest)
{
    char fname[FILENAME_MAX];
    int err;

    err = find_or_download_pdf(GRETL_GUIDE, get_manpref(), fname);

#ifdef G_OS_WIN32
    if (!err) {
	gretl_show_pdf(fname, dest);
    }
#elif defined(OS_OSX)
    if (!err) {
	gretl_show_pdf(fname, dest);
    }
#else /* Linux */
    if (!err) {
	gchar *tmp = NULL;

	if (strstr(viewpdf, "okular")) {
	    /* special case: option stuck onto fname */
	    tmp = g_strdup_printf("%s#%s", fname, dest);
	    gretl_show_pdf(tmp, NULL);
	} else {
	    if (strstr(viewpdf, "xpdf")) {
		tmp = g_strdup_printf("+%s", dest);
	    } else if (strstr(viewpdf, "evince")) {
		tmp = g_strdup_printf("--named-dest=%s", dest);
	    }
	    gretl_show_pdf(fname, tmp);
	}
	g_free(tmp);
    }
#endif
}

void display_gnuplot_help (void)
{
    char fname[FILENAME_MAX];
    int err;

    err = find_or_download_pdf(GNUPLOT_REF, 0, fname);

    if (!err) {
	gretl_show_pdf(fname, NULL);
    }
}

void display_x12a_help (void)
{
    char fname[FILENAME_MAX];
    int err;

    err = find_or_download_pdf(X12A_REF, 0, fname);

    if (!err) {
	gretl_show_pdf(fname, NULL);
    }
}
