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

/* session.c for gretl */

#include "gretl.h"
#include "gui_utils.h"
#include "session.h"
#include "selector.h"
#include "ssheet.h"
#include "plotspec.h"
#include "gpt_control.h"
#include "guiprint.h"
#include "gui_recode.h"
#include "model_table.h"
#include "graph_page.h"
#include "textbuf.h"
#include "cmdstack.h"
#include "filelists.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "menustate.h"
#include "toolbar.h"
#include "winstack.h"
#include "fncall.h"
#include "lib_private.h"

#include "var.h"
#include "varprint.h"
#include "objstack.h"
#include "system.h"
#include "gretl_xml.h"
#include "gretl_func.h"
#include "matrix_extra.h"
#include "uservar.h"
#include "gretl_zip.h"
#include "dbread.h"

#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>

#ifdef _WIN32
# include <windows.h>
# include "gretlwin32.h"
#endif

#include "../pixmaps/model.xpm"
#include "../pixmaps/boxplot.xpm"
#include "../pixmaps/gnuplot.xpm"
#include "../pixmaps/xfm_sc.xpm"
#include "../pixmaps/xfm_info.xpm"
#include "../pixmaps/xfm_text.xpm"
#include "../pixmaps/rhohat.xpm"
#include "../pixmaps/summary.xpm"
#include "../pixmaps/model_table.xpm"
#include "../pixmaps/graph_page.xpm"
#include "../pixmaps/matrix.xpm"
#include "../pixmaps/bundle.xpm"

#define SESSION_DEBUG 0
#define GRAPH_DEBUG 1
#define SESSION_BUNDLE 1

#define SHOWNAMELEN 15
#define ICONVIEW_MIN_COLS 4

#if GTK_MAJOR_VERSION < 3
# define REMEDY_LABELS
#endif

typedef struct SESSION_ SESSION;
typedef struct SESSION_TEXT_ SESSION_TEXT;
typedef struct SESSION_MODEL_ SESSION_MODEL;
typedef struct SESSION_GRAPH_ SESSION_GRAPH;
typedef struct gui_obj_ gui_obj;

enum {
    SESSION_CHANGED    = 1 << 0,
    SESSION_SAVED      = 1 << 1,
    SESSION_OPEN       = 1 << 2
};

struct SESSION_ {
    char name[MAXLEN];
    char dirname[MAXLEN];
    int status;
    int show_notes;
    int nmodels;
    int ngraphs;
    int ntexts;
    SESSION_GRAPH **graphs;
    SESSION_MODEL **models;
    SESSION_TEXT **texts;
    char *notes;
};

struct SESSION_TEXT_ {
    char name[MAXSAVENAME];
    char *buf;
};

struct SESSION_MODEL_ {
    char name[MAXSAVENAME];
    void *ptr;
    GretlObjType type;
};

struct SESSION_GRAPH_ {
    char name[MAXSAVENAME];
    char fname[MAXLEN];
    GretlObjType type;
    int has_datafile;
};

struct gui_obj_ {
    gchar *name;
    gint sort;
    gpointer data;
    GtkWidget *icon;
    GtkWidget *label;
    gint row, col;
};

struct sample_info {
    char datafile[MAXLEN];
    int t1;
    int t2;
    char *mask;
    char *restriction;
    unsigned int seed;
    int resample_n;
};

enum {
    ICON_ADD_BATCH,
    ICON_ADD_SINGLE
} icon_add_modes;

static char *global_items[] = {
    N_("Save session"),
    N_("Arrange icons"),
    N_("Add matrix..."),
    N_("Windows"),
    N_("Close window")
};

static char *model_items[] = {
    N_("Display"),
    N_("Add to model table"),
    N_("Rename"),
    N_("Delete")
};
#define ADD_TO_MTAB_IDX 1 /* position of "Add to model table" */

static char *model_table_items[] = {
    N_("Display"),
    N_("Clear"),
    N_("Help")
};

static char *graph_page_items[] = {
    N_("Display"),
    N_("Save as TeX..."),
    N_("Clear"),
    N_("Help")
};

static char *generic_items[] = {
    N_("Display"),
    N_("Rename"),
    N_("Delete")
};

static char *graph_items[] = {
    N_("Display"),
    N_("Edit plot commands"),
    N_("Add to graph page"),
    N_("Rename"),
    N_("Delete"),
    N_("Copy")
};
#define ADD_TO_GPAGE_IDX 2 /* position of "Add to graph page" */
#define GRAPH_COPY_IDX 5 /* position of "Copy" */

static char *dataset_items[] = {
    N_("Edit"),
    N_("Export as CSV..."),
    N_("Copy as CSV...")
};

static char *scalars_items[] = {
    N_("Edit"),
    N_("Copy as CSV...")
};

static char *info_items[] = {
    N_("View")
};

static char *matrix_items[] = {
    N_("View"),
    N_("Edit"),
    N_("Properties"),
    N_("Copy as CSV..."),
    N_("Rename"),
    N_("Delete"),
    N_("Copy")
};

static char *bundle_items[] = {
    N_("View"),
    N_("Rename"),
    N_("Delete")
};

/* file-scope globals */

SESSION session; /* holds named models, graphs, etc. */

static char sessionfile[MAXLEN];

static GtkWidget *iconview;
static GtkWidget *icon_table;
static GtkWidget *global_popup;
static GtkWidget *model_popup;
static GtkWidget *model_table_popup;
static GtkWidget *generic_popup;
static GtkWidget *graph_popup;
static GtkWidget *graph_page_popup;
static GtkWidget *data_popup;
static GtkWidget *scalars_popup;
static GtkWidget *info_popup;
static GtkWidget *matrix_popup;
static GtkWidget *bundle_popup;
static GtkWidget *save_item;
static GtkWidget *mtab_add_item;
static GtkWidget *gpage_add_item;
static GtkWidget *graph_copy_item;

static GList *iconlist;
static gui_obj *active_object;
static gint iconview_width;
static gint iconview_cols;
static gint in_icon;

/* private functions */
static gui_obj *gui_object_new (gchar *name, int sort, gpointer data);
static gui_obj *session_add_icon (gpointer data, int sort, int mode);
static void session_build_popups (void);
static void global_popup_callback (GtkWidget *widget, gpointer data);
static void object_popup_callback (GtkWidget *widget, gpointer data);
static void data_popup_callback (GtkWidget *widget, gpointer data);
static void scalars_popup_callback (GtkWidget *widget, gpointer data);
static void info_popup_callback (GtkWidget *widget, gpointer data);
static void matrix_popup_callback (GtkWidget *widget, gpointer data);
static void bundle_popup_callback (GtkWidget *widget, gpointer data);
static void session_delete_icon (gui_obj *obj);
static void open_gui_graph (gui_obj *obj);
static gboolean session_view_click (GtkWidget *widget,
				    GdkEventButton *event,
				    gpointer data);
static int real_delete_model_from_session (SESSION_MODEL *model);
static void make_short_label_string (char *targ, const char *src);
static gui_obj *get_gui_obj_by_data (void *finddata);
static int gui_user_var_callback (const char *name, GretlType type,
				  int action);
static void auto_view_session (void);
static int display_session_model (SESSION_MODEL *sm);

static int session_graph_count;
static int session_bundle_count;
static int commands_recorded;

int session_is_modified (void)
{
    return session.status & SESSION_CHANGED;
}

int session_is_open (void)
{
    return session.status & SESSION_OPEN;
}

void mark_session_changed (void)
{
    iconview_menubar_state(TRUE);
    session.status &= ~SESSION_SAVED;
    session.status |= SESSION_CHANGED;
    if (save_item != NULL) {
	gtk_widget_set_sensitive(save_item, TRUE);
    }
    flip(mdata->ui, "/menubar/File/SessionFiles/SaveSession", TRUE);
    flip(mdata->ui, "/menubar/File/SessionFiles/SaveSessionAs", TRUE);
    if (*session.name != '\0') {
	set_main_window_title(session.name, TRUE);
    }
}

static void mark_session_saved (void)
{
    session.status &= ~SESSION_CHANGED;
    session.status |= SESSION_SAVED;
    if (save_item != NULL) {
	gtk_widget_set_sensitive(save_item, FALSE);
    }
    flip(mdata->ui, "/menubar/File/SessionFiles/SaveSession", FALSE);
    flip(mdata->ui, "/menubar/File/SessionFiles/SaveSessionAs",
	 session.status & SESSION_OPEN);
    set_main_window_title(session.name, FALSE);
}

void set_commands_recorded (void)
{
    commands_recorded = 1;
}

int get_commands_recorded (void)
{
    return commands_recorded;
}

int widget_is_iconview (GtkWidget *w)
{
    return w == iconview;
}

/* constructors and destructors for session data-objects */

static void free_session_text (SESSION_TEXT *text)
{
    free(text->buf);
    free(text);
}

static void free_session_model (SESSION_MODEL *mod)
{
#if SESSION_DEBUG
    fprintf(stderr, "free_session_model: ptr at %p\n", (void *) mod->ptr);
#endif
    gretl_object_remove_from_stack(mod->ptr, mod->type);
    free(mod);
}

static SESSION_MODEL *session_model_new (void *ptr, const char *name,
					 GretlObjType type)
{
    SESSION_MODEL *mod = mymalloc(sizeof *mod);

    if (mod != NULL) {
	mod->ptr = ptr;
	mod->type = type;
	gretl_stack_object_as(ptr, type, name);
	strcpy(mod->name, gretl_object_get_name(ptr, type));
	/* note: we don't add a reference here: that's handled by the
	   mechanism in objstack.c
	*/
    }

    return mod;
}

static int session_append_text (const char *tname, char *buf)
{
    SESSION_TEXT *text;
    SESSION_TEXT **texts;
    int nt = session.ntexts;

    text = mymalloc(sizeof *text);
    if (text == NULL) {
	return 1;
    }

    text->buf = buf;

    *text->name = '\0';
    if (tname != NULL) {
	strncat(text->name, tname, MAXSAVENAME - 1);
    }

    texts = myrealloc(session.texts, (nt + 1) * sizeof *texts);

    if (texts == NULL) {
	free_session_text(text);
	return 1;
    }

    session.texts = texts;
    session.texts[nt] = text;
    session.ntexts += 1;
    mark_session_changed();

    return 0;
}

static int session_append_model (SESSION_MODEL *mod)
{
    SESSION_MODEL **models;
    int nm = session.nmodels;

    models = myrealloc(session.models, (nm + 1) * sizeof *models);
    if (models == NULL) {
	free_session_model(mod);
	return 1;
    }

    session.models = models;
    session.models[nm] = mod;
    session.nmodels += 1;

#if SESSION_DEBUG
    fprintf(stderr, "session_append_model: nmodels now = %d\n",
	    session.nmodels);
#endif

    return 0;
}

static SESSION_GRAPH *session_append_graph (const char *grname,
					    const char *fname,
					    GretlObjType type)
{
    SESSION_GRAPH *graph;
    SESSION_GRAPH **graphs;
    int ng = session.ngraphs;

    graph = mymalloc(sizeof *graph);
    if (graph == NULL) {
	return NULL;
    }

    *graph->name = '\0';
    if (grname != NULL) {
	strncat(graph->name, grname, MAXSAVENAME - 1);
    }

    *graph->fname = '\0';
    if (fname != NULL) {
	strncat(graph->fname, fname, MAXLEN - 1);
    }

    graph->type = type;
    graph->has_datafile = 0;

    graphs = myrealloc(session.graphs, (ng + 1) * sizeof *graphs);

    if (graphs == NULL) {
	free(graph);
	return NULL;
    }

    session.graphs = graphs;
    session.graphs[ng] = graph;
    session.ngraphs++;     /* decremented if graph is deleted */
    session_graph_count++; /* not decremented on deletion */

    return graph;
}

const char *last_session_graph_name (void)
{
    int ng = session.ngraphs;

    if (ng > 0) {
	SESSION_GRAPH *graph = session.graphs[ng - 1];

	return graph->fname;
    } else {
	return NULL;
    }
}

static void session_switch_log_location (int code)
{
    gchar *dirpath;

    dirpath = gretl_make_dotpath(session.dirname);
    set_session_log(dirpath, code);
    g_free(dirpath);
}

char *session_graph_make_path (char *path, const char *fname)
{
    if (strstr(fname, session.dirname) != NULL) {
	/* should be OK */
	strcpy(path, fname);
    } else {
	const char *p = path_last_slash_const(fname);

	if (p != NULL) {
	    sprintf(path, "%s%s", session.dirname, p);
	} else {
	    sprintf(path, "%s%c%s", session.dirname, SLASH, fname);
	}
    }

    return path;
}

/* first arg should be a MAXLEN string */

static char *session_file_make_path (char *path,
				     const char *fname,
				     const char *sdir)
{
#if SESSION_DEBUG
    if (sdir != NULL) {
	fprintf(stderr, "session_file_make_path: fname '%s', sdir '%s'\n",
		fname, sdir);
    } else {
	fprintf(stderr, "session_file_make_path: fname '%s', session.dirname '%s'\n",
		fname, session.dirname);
    }
#endif

    if (g_path_is_absolute(fname)) {
	strcpy(path, fname);
    } else if (sdir != NULL) {
	gretl_build_path(path, gretl_dotdir(), sdir, fname, NULL);
    } else {
	gretl_build_path(path, gretl_dotdir(), session.dirname,
			 fname, NULL);
    }

#if SESSION_DEBUG
    fprintf(stderr, "session_file_make_path: outgoing path '%s'\n", path);
#endif

    return path;
}

#include "session_xml.c"

static void edit_session_notes (void)
{
    static GtkWidget *notes_window;

    if (notes_window == NULL) {
	windata_t *vwin;

	vwin = edit_buffer(&session.notes, 80, 400,
			   _("gretl: session notes"),
			   EDIT_NOTES);
	notes_window = vwin->main;
	g_signal_connect(G_OBJECT(vwin->main), "destroy",
			 G_CALLBACK(gtk_widget_destroyed),
			 &notes_window);
    } else {
	gtk_window_present(GTK_WINDOW(notes_window));
    }
}

int is_session_model (void *p)
{
    int i;

#if SESSION_DEBUG
    fprintf(stderr, "is_session_model: testing %p (nmodels = %d)\n",
	    p, session.nmodels);
#endif

    for (i=0; i<session.nmodels; i++) {
	if (p == session.models[i]->ptr) {
	    return 1;
	}
    }
    return 0;
}

static SESSION_MODEL *get_session_model_by_name (const char *name)
{
    int i;

    for (i=0; i<session.nmodels; i++) {
	if (!strcmp(name, session.models[i]->name)) {
	    return session.models[i];
	}
    }

    return NULL;
}

static SESSION_GRAPH *get_session_graph_by_name (const char *name)
{
    int i;

    for (i=0; i<session.ngraphs; i++) {
	if (!strcmp(name, session.graphs[i]->name)) {
	    return session.graphs[i];
	}
    }

    return NULL;
}

static SESSION_TEXT *get_session_text_by_name (const char *name)
{
    int i;

    for (i=0; i<session.ntexts; i++) {
	if (!strcmp(name, session.texts[i]->name)) {
	    return session.texts[i];
	}
    }

    return NULL;
}

int real_add_text_to_session (PRN *prn, int pos, const char *tname)
{
    SESSION_TEXT *text = get_session_text_by_name(tname);
    char *buf = gretl_print_get_chunk_at(prn, pos);
    int replace = 0;

    if (buf == NULL || string_is_blank(buf)) {
	warnbox_printf(_("%s: no text to save"), tname);
	return ADD_OBJECT_FAIL;
    }

    if (text != NULL) {
	free(text->buf);
	text->buf = buf;
	replace = 1;
    } else {
	if (session_append_text(tname, buf)) {
	    return ADD_OBJECT_FAIL;
	}
    }

    mark_session_changed();

    if (iconlist != NULL && !replace) {
	session_add_icon(text, GRETL_OBJ_TEXT, ICON_ADD_SINGLE);
    } else if (autoicon_on()) {
	auto_view_session();
    }

    return (replace)? ADD_OBJECT_REPLACE : ADD_OBJECT_OK;
}

static int ends_with_number (const char *s, int *k)
{
    int ret = 0;

    s = strrchr(s, '(');

    if (s != NULL) {
	int n = strspn(s + 1, "0123456789");

	if (n > 0 && *(s + n + 1) == ')' && *(s + n + 2) == '\0') {
	    sscanf(s + 1, "%d", k);
	    ret = 1;
	}
    }

    return ret;
}

static void ensure_unique_text_name (char *tname)
{
    char tmp[MAXSAVENAME];
    const char *p, *oldname;
    int i, n, id, idmax = 0;

    for (i=0; i<session.ntexts; i++) {
	id = 0;
	oldname = session.texts[i]->name;
	if (!strcmp(tname, oldname)) {
	    id = 1;
	} else if (ends_with_number(oldname, &id)) {
	    p = strrchr(oldname, '(');
	    *tmp = '\0';
	    strncat(tmp, oldname, p - oldname);
	    if (strcmp(tmp, tname)) {
		id = 0;
	    }
	}
	if (id > idmax) {
	    idmax = id;
	}
    }

    if (idmax > 0) {
	char num[16];

	sprintf(num, "(%d)", ++idmax);
	n = MAXSAVENAME - strlen(num) - 1;
	*tmp = '\0';
	strncat(tmp, tname, n);
	strcat(tmp, num);
	strcpy(tname, tmp);
    }
}

void save_output_as_text_icon (windata_t *vwin)
{
    const gchar *title = NULL;
    char tname[MAXSAVENAME];
    gchar *buf;
    int err;

    buf = textview_get_text(vwin->text);
    if (buf == NULL) {
	errbox("Couldn't retrieve buffer");
	return;
    }

    if (vwin->main != NULL) {
	title = gtk_window_get_title(GTK_WINDOW(vwin->main));
    }

    if (title != NULL && !strncmp(title, "gretl: ", 7)) {
	*tname = '\0';
	strncat(tname, title + 7, MAXSAVENAME - 1);
    } else {
	strcpy(tname, "text");
    }

    ensure_unique_text_name(tname);
    err = session_append_text(tname, buf);

    if (err) {
	return;
    } else if (iconlist != NULL) {
	SESSION_TEXT * text = get_session_text_by_name(tname);

	session_add_icon(text, GRETL_OBJ_TEXT, ICON_ADD_SINGLE);
    } else if (autoicon_on()) {
	auto_view_session();
    }
}

static int add_model_to_session (void *ptr, const char *name,
				 GretlObjType type,
				 SESSION_MODEL **psm)
{
    SESSION_MODEL *model;
    int err = 0;

#if SESSION_DEBUG
    fprintf(stderr, "add_model_to_session: doing session_model_new\n"
	    " with ptr = %p\n", ptr);
#endif

    model = session_model_new(ptr, name, type);

    if (model == NULL || session_append_model(model)) {
	err = E_ALLOC;
    } else if (iconlist != NULL) {
	session_add_icon(model, type, ICON_ADD_SINGLE);
    } else if (autoicon_on()) {
	auto_view_session();
    }

    if (!err && psm != NULL) {
	*psm = model;
    }

    return err;
}

static int replace_session_graph (SESSION_GRAPH *graph,
				  const char *fname,
				  GretlObjType type)
{
    if (fname != NULL && strcmp(graph->fname, fname)) {
	gretl_remove(graph->fname);
	strcpy(graph->fname, fname);
    }

    graph->type = type;
    mark_session_changed();

    return ADD_OBJECT_REPLACE;
}

static int real_add_graph_to_session (const char *fname,
				      const char *grname,
				      GretlObjType type,
				      SESSION_GRAPH **pgraph)
{
    SESSION_GRAPH *graph = get_session_graph_by_name(grname);
    int ret = ADD_OBJECT_OK;

    if (graph != NULL) {
	ret = replace_session_graph(graph, grname, type);
    } else {
	graph = session_append_graph(grname, fname, type);
	if (graph == NULL) {
	    ret = ADD_OBJECT_FAIL;
	}
    }

    if (ret != ADD_OBJECT_FAIL) {
	if (pgraph != NULL) {
	    *pgraph = graph;
	}
	mark_session_changed();
	if (iconlist != NULL) {
	    session_add_icon(graph, type, ICON_ADD_SINGLE);
	    if (autoicon_on()) {
		gtk_window_present(GTK_WINDOW(iconview));
	    }
	} else if (autoicon_on()) {
	    auto_view_session();
	}
    }

    return ret;
}

const char *get_session_dirname (void)
{
    return session.dirname;
}

static int session_dir_ok (void)
{
    int ret = 1;

    if (*session.dirname == '\0') {
	pid_t p = getpid();

	errno = 0;

	sprintf(session.dirname, ".gretl-%d", (int) p);
	if (gretl_chdir(gretl_dotdir())) {
	    perror("moving to user directory");
	    ret = 0;
	} else if (gretl_mkdir(session.dirname)) {
	    ret = 0;
	}
    }

    return ret;
}

static void make_graph_name (char *shortname, char *graphname)
{
    int i, n = session_graph_count + 1;

    for ( ; ; n++) {
	int unique = 1;

	sprintf(shortname, "graph.%d", n);
	sprintf(graphname, "%s %d", _("Graph"), n);

	for (i=0; i<session.ngraphs; i++) {
	    if (!strcmp(shortname, session.graphs[i]->fname) ||
		!strcmp(graphname, session.graphs[i]->name)) {
		unique = 0;
		break;
	    }
	}

	if (unique) {
	    break;
	}
    }
}

static int maybe_move_plot_datafile (const char *orig,
				     const char *revised,
				     int *has_datafile)
{
    gchar *tmp1 = g_strdup_printf("%s.dat", orig);
    int err = 0;

    if (gretl_test_fopen(tmp1, "r") == 0) {
	gchar *tmp2 = g_strdup_printf("%s.dat", revised);

	err = gretl_rename(tmp1, tmp2);
	g_free(tmp2);
	*has_datafile = 1;
    }

    g_free(tmp1);

    return err;
}

/* Callback for "add to session as icon" on a graph displayed
   as PNG -- see gpt_control.c. Note that there is code in
   gpt_control designed to ensure that this option is not
   available for a graph that has already been saved in
   this way.
*/

int gui_add_graph_to_session (char *fname, char *fullname, int type)
{
    char shortname[MAXSAVENAME];
    char graphname[MAXSAVENAME];
    int has_datafile = 0;
    int err = 0;

    if (type != GRETL_OBJ_PLOT && (dataset == NULL || dataset->Z == NULL)) {
	/* we may be called via the "stats calculator" when
	   there's no dataset yet */
	err = open_nulldata(dataset, 0, 10, OPT_NONE, NULL);
	if (err) {
	    gui_errmsg(err);
	    return 1;
	}
	register_data(NULLDATA_STARTED);
    }

    errno = 0;

    if (!session_dir_ok()) {
	errbox(_("Failed to copy graph file"));
	return 1;
    }

    err = gretl_chdir(gretl_dotdir());
    if (err) {
	gui_errmsg(err);
    } else {
	make_graph_name(shortname, graphname);
	session_file_make_path(fullname, shortname, NULL);

	/* copy temporary plot file to session directory */
	err = copyfile(fname, fullname);
	if (!err) {
	    /* remove the original and transcribe the new
	       name to @fname */
	    maybe_move_plot_datafile(fname, fullname, &has_datafile);
	    gretl_remove(fname);
	    strcpy(fname, fullname);
	}
    }

    if (!err) {
	SESSION_GRAPH *grf = NULL;

	err = real_add_graph_to_session(shortname, graphname, type, &grf);
	if (err == ADD_OBJECT_FAIL) {
	    err = 1;
	} else if (has_datafile) {
	    grf->has_datafile = 1;
	}
    }

    return err;
}

static void make_graph_filename (char *shortname)
{
    int i, n = session_graph_count + 1;

    for ( ; ; n++) {
	int unique = 1;

	sprintf(shortname, "graph.%d", n);

	for (i=0; i<session.ngraphs; i++) {
	    if (!strcmp(shortname, session.graphs[i]->fname)) {
		unique = 0;
		break;
	    }
	}

	if (unique) {
	    break;
	}
    }
}

/* Callback for a command on the pattern "gname <- plotspec".
   Note @gname may contain non-ASCII characters (which will be
   in UTF-8). So we'd best not use @gname to construct a filename
   for the graph in case we're on a non-UTF-8 system.
*/

int cli_add_graph_to_session (const char *fname, const char *gname,
			      GretlObjType type, int display)
{
    SESSION_GRAPH *graph = NULL;
    char shortname[MAXSAVENAME];
    char grpath[MAXLEN];
    int replace = 0;
    int ret = 0;

    errno = 0;

    if (!session_dir_ok()) {
	errbox(_("Failed to copy graph file"));
	return ADD_OBJECT_FAIL;
    }

    gretl_chdir(gretl_dotdir());

    graph = get_session_graph_by_name(gname);

    if (graph != NULL) {
	/* replacing */
	session_file_make_path(grpath, graph->fname, NULL);
	replace = 1;
    } else {
	/* ensure unique filename */
	make_graph_filename(shortname);
	session_file_make_path(grpath, shortname, NULL);
    }

    if (copyfile(fname, grpath)) {
	errbox(_("Failed to copy graph file"));
	return ADD_OBJECT_FAIL;
    }

    /* we copied the plot commands file, @fname, into the
       session directory; now delete the original */
    gretl_remove(fname);

    if (replace) {
	ret = replace_session_graph(graph, NULL, type);
    } else {
	ret = real_add_graph_to_session(shortname, gname, type, &graph);
    }

    if (ret == ADD_OBJECT_REPLACE) {
	/* don't keep a stale plot window open */
	GtkWidget *pwin = get_window_for_plot(graph);

	if (pwin != NULL) {
	    gtk_widget_destroy(pwin);
	}
    }

    if (ret != ADD_OBJECT_FAIL && display) {
	display_session_graph_by_data(graph);
    }

    return ret;
}

void *get_session_object_by_name (const char *name, GretlObjType *type)
{
    int i;

    for (i=0; i<session.nmodels; i++) {
	if (!strcmp(name, session.models[i]->name)) {
	    *type = session.models[i]->type;
	    return session.models[i]->ptr;
	}
    }

    for (i=0; i<session.ngraphs; i++) {
	if (!strcmp(name, session.graphs[i]->name)) {
	    *type = GRETL_OBJ_GRAPH;
	    return session.graphs[i];
	}
    }

    for (i=0; i<session.ntexts; i++) {
	if (!strcmp(name, session.texts[i]->name)) {
	    *type = GRETL_OBJ_TEXT;
	    return session.texts[i];
	}
    }

    return NULL;
}

static SESSION_MODEL *get_session_model_by_data (const void *ptr)
{
    int i;

    for (i=0; i<session.nmodels; i++) {
	if (ptr == session.models[i]->ptr) {
	    return session.models[i];
	}
    }

    return NULL;
}

static void delete_icon_for_data (void *data)
{
    gui_obj *gobj = get_gui_obj_by_data(data);

    if (gobj != NULL) {
	/* icon is shown: delete it */
	session_delete_icon(gobj);
    }
}

static int show_model_in_window (void *ptr, GretlObjType type)
{
    char *name = NULL;
    PRN *prn;

    if (type != GRETL_OBJ_EQN &&
	type != GRETL_OBJ_VAR &&
	type != GRETL_OBJ_SYS) {
	return 1;
    }

    if (bufopen(&prn)) {
	return 1;
    }

    gretl_object_compose_unique_name(ptr, type);
    name = gretl_object_get_name(ptr, type);

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) ptr;

	printmodel(pmod, dataset, OPT_NONE, prn);
	view_model(prn, pmod, name);
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;

	gretl_VAR_print(var, dataset, OPT_NONE, prn);
	view_buffer(prn, 78, 450, name, var->ci, var);
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) ptr;

	gretl_system_print(sys, dataset, OPT_NONE, prn);
	view_buffer(prn, 78, 450, name, SYSTEM, sys);
    }

    return 0;
}

/* Callback (indirectly) from libgretl for a model created via
   script command. When this function is called, the model
   in question has already been stacked; it is just a matter
   of syncing the GUI session with the model stack state.
*/

int add_model_to_session_callback (void *ptr, GretlObjType type,
				   gretlopt opt)
{
    SESSION_MODEL *model = NULL;
    char *name = gretl_object_get_name(ptr, type);
    int err = 0;

    if (type == GRETL_OBJ_SYS && (opt & OPT_W)) {
	/* we got the --window callback from "estimate" */
	return show_model_in_window(ptr, type);
    } else if (name == NULL || *name == '\0') {
	/* we just got the --window callback */
	return show_model_in_window(ptr, type);
    }

    /* are we replacing a session model's content? */
    model = get_session_model_by_name(name);

    if (model != NULL) {
	windata_t *vwin = get_viewer_for_data(model->ptr);

	if (vwin != NULL) {
	    /* current model window is "orphaned" */
	    gretl_viewer_destroy(vwin);
	}
	model->ptr = ptr;
	model->type = type;
	mark_session_changed();
    } else {
	err = add_model_to_session(ptr, name, type, &model);
	if (!err) {
	    mark_session_changed();
	}
    }

    if (!err && (opt & OPT_W)) {
	/* should open a window for this session model */
	display_session_model(model);
    }

    return err;
}

static int model_type_from_vwin (windata_t *vwin)
{
    if (vwin->role == SYSTEM) {
	return GRETL_OBJ_SYS;
    } else if (vwin->role == VAR || vwin->role == VECM) {
	return GRETL_OBJ_VAR;
    } else {
	return GRETL_OBJ_EQN;
    }
}

static int close_on_add (GtkAction *action)
{
    if (action == NULL) {
	return 1;
    } else {
	return (strstr(gtk_action_get_name(action), "Close") != NULL);
    }
}

void model_add_as_icon (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    void *ptr = vwin->data;
    const char *name;
    int type;
    int err = 0;

    if (ptr == NULL) {
	return;
    }

    gretl_error_clear();

#if SESSION_DEBUG
    fprintf(stderr, "model_add_as_icon: ptr = %p\n", ptr);
#endif

    if (get_session_model_by_data(ptr)) {
	/* "can't happen" */
	if (close_on_add(action)) {
	    gretl_viewer_destroy(vwin);
	} else {
	    infobox(_("Model is already saved"));
	}
	return;
    }

    type = model_type_from_vwin(vwin);
    name = gretl_object_get_name(ptr, type);

    if (name != NULL && *name == '\0') {
	/* since we didn't get a match by data above, we
	   need to ensure that this model is given a
	   unique name */
	gretl_object_compose_unique_name(ptr, type);
	name = gretl_object_get_name(ptr, type);
    }

    if (name == NULL) {
	nomem();
	return;
    }

    err = add_model_to_session(ptr, name, type, NULL);

    if (!err) {
	mark_session_changed();
	if (close_on_add(action)) {
	    gretl_viewer_destroy(vwin);
	} else {
	    set_model_save_state(vwin, FALSE);
	}
    }
}

/* Called (via gui_utils.c) from toolbar.c, to implement
   saving displayed bundle as icon; handles VIEW_BUNDLE
   and also VIEW_DBNOMICS
*/

void bundle_add_as_icon (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    gretl_bundle *bundle = vwin->data;
    char vname[VNAMELEN];
    gchar *blurb;
    int resp;

    sprintf(vname, "bundle%d", ++session_bundle_count);
    blurb = g_strdup_printf("Save bundle\nName (max. %d characters):",
			    VNAMELEN - 1);
    resp = object_name_entry_dialog(vname, GRETL_TYPE_BUNDLE,
				    blurb, NULL, vwin->main);
    g_free(blurb);

    if (!canceled(resp)) {
	int err = user_var_add(vname, GRETL_TYPE_BUNDLE, bundle);

	if (err) {
	    gui_errmsg(err);
	} else {
	    mark_session_changed();
	    vwin_action_set_sensitive(vwin, "SaveAsIcon", FALSE);
	}
    }
}

static void
session_name_from_session_file (char *sname, const char *fname)
{
    const char *p = path_last_slash_const(fname);
    char *q;

    if (p != NULL) {
	strcpy(sname, p + 1);
    } else {
	strcpy(sname, fname);
    }

    q = strstr(sname, ".gretl");
    if (q != NULL && strlen(q) == 6) {
	*q = '\0';
    }

#if SESSION_DEBUG
    fprintf(stderr, "session_name_from_session_file: %s -> %s\n",
	    fname, sname);
#endif
}

/* remedial action in case of mis-coded filename */

static int get_session_dataname (char *fname, const char *sdir)
{
    const gchar *dname;
    GDir *dir;
    FILE *fp;
    gchar *tmp;
    int n, err = E_FOPEN;

    if (sdir != NULL) {
	tmp = g_strdup(sdir);
    } else {
	tmp = g_strdup(session.dirname);
    }
    n = strlen(tmp);

    if (tmp[n-1] == '/' || tmp[n-1] == '\\') {
	tmp[n-1] = '\0';
    }

    dir = gretl_opendir(tmp);

    if (dir != NULL) {
	while ((dname = g_dir_read_name(dir)) != NULL) {
	    if (has_suffix(dname, ".gdt")) {
		session_file_make_path(fname, dname, sdir);
		fp = gretl_fopen(fname, "r");
		if (fp != NULL) {
		    fclose(fp);
		    err = 0;
		}
		break;
	    }
	}
	g_dir_close(dir);
    }

    g_free(tmp);

    return err;
}

static void sinfo_init (struct sample_info *sinfo)
{
    strcpy(sinfo->datafile, "data.gdt");
    sinfo->t1 = 0;
    sinfo->t2 = 0;
    sinfo->mask = NULL;
    sinfo->restriction = NULL;
    sinfo->seed = 0;
    sinfo->resample_n = 0;
}

static void sinfo_free_data (struct sample_info *sinfo)
{
    if (sinfo->mask != NULL) {
	free(sinfo->mask);
    }
    if (sinfo->restriction != NULL) {
	free(sinfo->restriction);
    }
}

static int test_session_dirname (const char *zdirname)
{
    char test[2*MAXLEN];

    g_return_val_if_fail(zdirname != NULL, 1);

    sprintf(test, "%s%csession.xml", zdirname, SLASH);

    if (gretl_test_fopen(test, "r") != 0) {
	return E_FOPEN;
    } else {
	return 0;
    }
}

static char *maybe_absolutize_tryfile (void)
{
    char *fname = get_tryfile();

    if (!g_path_is_absolute(fname)) {
	gchar *cwd = g_get_current_dir();

	if (cwd != NULL) {
	    gchar *tmp = g_build_filename(cwd, fname, NULL);

	    set_tryfile(tmp);
	    g_free(tmp);
	    g_free(cwd);
	}
    }

    return get_tryfile();
}

/* note: the name of the file to be opened is in the gretl.c var
   'tryfile' */

static gboolean real_open_session (gretl_bundle **pb)
{
    struct sample_info sinfo;
    char xmlname[MAXLEN]; /* path to master session XML file */
    char gdtname[MAXLEN]; /* path to session data file */
    char fname[MAXLEN];   /* multi-purpose temp variable */
    char *tryname = get_tryfile();
    gchar *zdirname = NULL;
    DATASET *sdset = NULL;
    FILE *fp;
    int as_bundle = 0;
    int nodata = 0;
    int err = 0;

    sinfo_init(&sinfo);

    fp = gretl_fopen(tryname, "r");
    if (fp != NULL) {
	fclose(fp);
    } else {
	file_read_errbox(tryname);
	delete_from_filelist(FILE_LIST_SESSION, tryname);
	return FALSE;
    }

    as_bundle = (pb != NULL);

    if (as_bundle) {
	sdset = datainfo_new();
    } else {
	/* close existing session, if any, and initialize */
	close_session(OPT_NONE);
	sdset = dataset;
    }

#if SESSION_DEBUG
    fprintf(stderr, "\nReading session file %s\n", tryname);
#endif

    /* we're about to change directory: if tryfile is not
       an absolute path we'll lose track of it
    */
    tryname = maybe_absolutize_tryfile();
    gretl_chdir(gretl_dotdir());
    err = gretl_unzip_session_file(tryname, &zdirname);

    if (err) {
	gui_errmsg(err);
	goto bailout;
    }

    if (!as_bundle) {
	session_name_from_session_file(session.name, tryname);
    }

    err = test_session_dirname(zdirname);

    if (err) {
	fprintf(stderr, "Failed on test_session_dirname\n");
	file_read_errbox("session.xml");
	goto bailout;
    } else {
	fprintf(stderr, "zdirname: '%s', OK\n", zdirname);
    }

    session_file_make_path(xmlname, "session.xml", zdirname);

    /* try getting the name of the session data file first */
    err = get_session_datafile_name(xmlname, &sinfo, &nodata);

    if (err) {
	fprintf(stderr, "get_session_datafile_name: err = %d\n", err);
	file_read_errbox("session.xml");
	goto bailout;
    }

    if (!nodata) {
	/* construct path to session data file */
	session_file_make_path(datafile, sinfo.datafile, zdirname);
	fp = gretl_fopen(datafile, "r");

	if (fp != NULL) {
	    /* OK, write good name into gdtname */
#if SESSION_DEBUG
	    fprintf(stderr, "got datafile name '%s'\n", datafile);
#endif
	    strcpy(gdtname, datafile);
	    fclose(fp);
	} else {
	    /* try remedial action, transform filename? */
	    fprintf(stderr, "'%s' : not found, trying to fix\n", datafile);
	    err = get_session_dataname(gdtname, zdirname);
	}

	if (!err) {
	    err = gretl_read_gdt(gdtname, sdset, OPT_B, NULL);
	}

	if (err) {
	    /* FIXME more explicit error message? */
	    file_read_errbox(sinfo.datafile);
	    goto bailout;
	} else {
#if SESSION_DEBUG
	    fprintf(stderr, "Opened session datafile '%s'\n", gdtname);
#endif
	    if (!as_bundle) {
		data_status = USER_DATA;
	    }
	}
    }

    if (as_bundle) {
	*pb = session_xml_to_bundle(xmlname, zdirname, sdset, &err);
	fprintf(stderr, "session_xml_to_bundle: err = %d\n", err);
	goto bailout;
    }

    /* having opened the data file (or not, if there's none), get the
       rest of the info from session.xml
    */
    strcpy(session.dirname, zdirname);
    err = read_session_xml(xmlname, &sinfo);
    if (err) {
	fprintf(stderr, "Failed on read_session_xml: err = %d\n", err);
	file_read_errbox("session.xml");
	goto bailout;
    }

    err = deserialize_user_vars(session.dirname);

    session_file_make_path(fname, "functions.xml", NULL);
    err = maybe_read_functions_file(fname);

    session_file_make_path(fname, "settings.inp", NULL);
    err = maybe_read_settings_file(fname);

    if (!nodata) {
	if (sinfo.resample_n > 0) {
	    err = dataset_resample(dataset, sinfo.resample_n, sinfo.seed);
	} else if (sinfo.mask != NULL) {
	    err = restrict_sample_from_mask(sinfo.mask, dataset, OPT_NONE);
	    if (!err) {
		dataset->restriction = sinfo.restriction;
		sinfo.restriction = NULL;
	    }
	}
	if (err) {
	    errbox(_("Couldn't set sample"));
	    goto bailout;
	} else {
	    dataset->t1 = sinfo.t1;
	    dataset->t2 = sinfo.t2;
	    register_data(OPENED_VIA_SESSION);
	    if (sinfo.mask != NULL) {
		set_sample_label(dataset);
	    }
	    sinfo_free_data(&sinfo);
	}
    }

    set_main_window_title(session.name, FALSE);

 bailout:

    g_free(zdirname);

    if (err) {
	delete_from_filelist(FILE_LIST_SESSION, tryname);
    } else if (!as_bundle) {
	strcpy(sessionfile, tryname);
	mkfilelist(FILE_LIST_SESSION, sessionfile, 0);

	session.status = SESSION_OPEN;
	/* sync gui with session */
	session_menu_state(TRUE);
	view_session();
	mark_session_saved();
	session_switch_log_location(LOG_OPEN);
	if (session.show_notes) {
	    edit_session_notes();
	}
    }

    if (as_bundle) {
	destroy_dataset(sdset);
    }

#if SESSION_DEBUG
    fprintf(stderr, "do_open_session: returning %d\n", !err);
#endif

    return !err;
}

gboolean do_open_session (void)
{
    return real_open_session(NULL);
}

gretl_bundle *open_session_as_bundle (void)
{
    gretl_bundle *ret = NULL;

    real_open_session(&ret);
    return ret;
}

void verify_clear_data (void)
{
    if (dataset_locked()) {
	return;
    }

    if (yes_no_dialog("gretl",
		      _("Clearing the data set will end\n"
			"your current session.  Continue?"),
		      NULL) != GRETL_YES) {
	return;
    }

    close_session(OPT_NONE); /* FIXME opt? */
}

static int remove_session_dir (void)
{
    gchar *fullpath;
    int err;

    fullpath = gretl_make_dotpath(session.dirname);
    err = gretl_chdir(gretl_dotdir());
    if (!err) {
	err = gretl_deltree(fullpath);
    }
    g_free(fullpath);

    return err;
}

void session_init (void)
{
    session.models = NULL;
    session.graphs = NULL;
    session.texts = NULL;
    session.notes = NULL;

    session.status = 0; /* note: neither changed nor saved */
    session.nmodels = 0;
    session.ngraphs = 0;
    session.ntexts = 0;

    *session.name = '\0';
    *session.dirname = '\0';

    set_user_var_callback(gui_user_var_callback);
}

void free_session (int on_exit)
{
    int i;

    if (session.models) {
	for (i=0; i<session.nmodels; i++) {
	    free_session_model(session.models[i]);
	}
	free(session.models);
	session.models = NULL;
    }
    session.nmodels = 0;

    if (session.graphs) {
	for (i=0; i<session.ngraphs; i++) {
	    free(session.graphs[i]);
	}
	free(session.graphs);
	session.graphs = NULL;
    }
    session.ngraphs = 0;

    if (session.texts) {
	for (i=0; i<session.ntexts; i++) {
	    free_session_text(session.texts[i]);
	}
	free(session.texts);
	session.texts = NULL;
    }
    session.ntexts = 0;

    if (session.notes) {
	free(session.notes);
	session.notes = NULL;
    }

    *session.name = '\0';

    if (*session.dirname != '\0') {
	if (on_exit) {
	    set_session_log(NULL, LOG_CLOSE);
	}
	remove_session_dir();
	*session.dirname = '\0';
    }
}

int highest_numbered_variable_in_session (void)
{
    GretlObjType type;
    void *ptr;
    int i, mvm, vmax = 0;

    if (session.models == NULL) {
	return 0;
    }

    for (i=0; i<session.nmodels; i++) {
	ptr = session.models[i]->ptr;
	if (ptr == NULL) {
	    continue;
	}
	type = session.models[i]->type;
	if (type == GRETL_OBJ_EQN) {
	    mvm = highest_numbered_var_in_model((MODEL *) ptr, dataset);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	} else if (type == GRETL_OBJ_VAR) {
	    mvm = gretl_VAR_get_highest_variable((GRETL_VAR *) ptr);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	} else if (type == GRETL_OBJ_SYS) {
	    mvm = highest_numbered_var_in_system((equation_system *) ptr,
						 dataset);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	}
    }

    return vmax;
}

GList *session_model_list (void)
{
    GList *list = NULL;

    if (session.models != NULL) {
	GretlObjType type;
	void *ptr;
	int i;

	for (i=0; i<session.nmodels; i++) {
	    ptr = session.models[i]->ptr;
	    type = session.models[i]->type;
	    if (ptr != NULL && type == GRETL_OBJ_EQN) {
		list = g_list_append(list, ptr);
	    }
	}
    }

    return list;
}

int session_file_is_open (void)
{
    return (*sessionfile != '\0');
}

void gui_clear_dataset (void)
{
    *datafile = 0;

    if (dataset->Z != NULL) {
	free_Z(dataset);
	dataset->Z = NULL;
    }

    clear_datainfo(dataset, CLEAR_FULL);

    clear_varlist(mdata->listbox);
    clear_sample_label();
    set_db_name(NULL, 0, NULL);

    data_status = 0;
    orig_vars = 0;
    dataset_menubar_state(FALSE);
}

static void session_clear_data (DATASET *pdinfo)
{
    gui_restore_sample(pdinfo);
    gui_clear_dataset();

    /* clear protected model */
    clear_model(model);

    free_command_stack();
    set_model_count(0);
    lib_cmd_destroy_context();
}

void close_session (gretlopt opt)
{
    int preserve = (opt & OPT_P)? 1 : 0;
    int logcode = LOG_NULL;
    int iview = 0;

#if SESSION_DEBUG
    fprintf(stderr, "close_session: starting cleanup\n");
#endif

    if (dataset != NULL && dataset->v > 0) {
	logcode = LOG_CLOSE;
	session_clear_data(dataset);
    }

    free_session(0);

    clear_model_table(1, NULL);
    clear_graph_page(1);

    session_menu_state(FALSE);
    *scriptfile = '\0';
    *sessionfile = '\0';

    if (iconview != NULL) {
	iview = 1;
	gtk_widget_destroy(iconview);
    }

    session.status = 0; /* not saved, changed or open */
    session.show_notes = 0;
    commands_recorded = 0;

    close_session_windows(opt);
    selector_cleanup();
    function_call_cleanup();
    edit_dialog_special_get_text(NULL);

    if (preserve) {
	/* preserve non-dataset items */
	libgretl_session_cleanup(SESSION_CLEAR_DATASET);
    } else {
	libgretl_session_cleanup(SESSION_CLEAR_ALL);
    }

    session_graph_count = 0;
    session_bundle_count = 0;
    reset_plot_count();
    reset_collection_count();

    set_session_log(NULL, logcode);

    if (iview && have_session_objects()) {
	dataset_menubar_state(FALSE);
	view_session();
    }
}

static void relpath_from_fname (char *path, const char *fname)
{
    const char *p;

    strcpy(path, ".");

    p = IS_SLASH(*fname) ? fname + 1 : fname;
    p = path_last_slash_const(p);

    if (p != NULL) {
	strcat(path, p + 1);
    } else {
	strcat(path, fname);
    }
}

/* dump the current dataset into the session dir */

static int real_save_session_dataset (const char *dname)
{
    char tmpname[MAXLEN];
    char *mask = NULL;
    char *restr = NULL;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    int write_err = 0;
    int err = 0;

    /* we need to retrieve and save the full version of the dataset */

    if (complex_subsampled()) {
	mask = copy_dataset_submask(dataset, &err);
	if (!err && dataset_is_resampled(dataset)) {
	    /* can't happen? */
	    mask = NULL;
	}
	if (dataset->restriction != NULL) {
	    restr = gretl_strdup(dataset->restriction);
	}
    }

    if (!err) {
	err = restore_full_sample(dataset, NULL);
    }

    if (!err) {
	session_file_make_path(tmpname, dname, NULL);
	write_err = gretl_write_gdt(tmpname, NULL, dataset,
				    OPT_NONE, 1);
    }

    if (mask != NULL) {
	/* reset the prior subsample */
	if (!err) {
	    err = restrict_sample_from_mask(mask, dataset,
					    OPT_NONE);
	}
	free(mask);
    }
    if (restr != NULL) {
	dataset->restriction = restr;
    }

    dataset->t1 = save_t1;
    dataset->t2 = save_t2;

    if (!err) {
	err = write_err;
    }

    if (!err) {
	/* flag the fact that the data are saved */
	data_status &= ~MODIFIED_DATA;
    }

    return err;
}

static const char *unpath (const char *fname)
{
    int i, n = strlen(fname);

    for (i=n-1; i>=0; i--) {
	if (fname[i] == '/') {
	    return fname + i + 1;
	}
#ifdef G_OS_WIN32
	if (fname[i] == '\\') {
	    return fname + i + 1;
	}
#endif
    }

    return fname;
}

/* called in the context of a re-opened session */

int save_session_dataset (void)
{
    int err = real_save_session_dataset(unpath(datafile));

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static void make_session_dataname (char *datname)
{
    int changed = 0;

    if (*datafile != '\0') {
	char *p;

	strcpy(datname, unpath(datafile));
	p = strrchr(datname, '.');
	if (p == NULL) {
	    strcat(datname, ".gdt");
	    changed = 1;
	} else if (strcmp(p, ".gdt")) {
	    strcpy(p, ".gdt");
	    changed = 1;
	}
    } else {
	strcpy(datname, "data.gdt");
	changed = 1;
    }

    if (changed) {
	GtkWidget *dlabel;

	dlabel = g_object_get_data(G_OBJECT(mdata->main), "dlabel");
	if (dlabel != NULL) {
	    gtk_label_set_text(GTK_LABEL(dlabel), datname);
	}
    }
}

#define SAVE_DEBUG 0

int save_session (char *fname)
{
    char *dirbak = NULL;
    char datname[MAXLEN];
    char dirname[MAXLEN];
    int len, err = 0;
    int log_code = LOG_SAVE;

    if (fname == NULL) {
	/* re-saving session 'as is' */
	fname = sessionfile;
    }

#if SAVE_DEBUG
    if (fname == sessionfile) {
	fprintf(stderr, "save_session:\n sessionfile='%s'\n",
		sessionfile);
    } else {
	fprintf(stderr, "save_session as:\n current session='%s'\n"
		" save filename='%s'\n", sessionfile, fname);
    }
#endif

    if (!session_dir_ok()) {
	errbox("Couldn't make session directory");
	return 1;
    }

    /* organize directory and file names */
    relpath_from_fname(dirname, fname);
    len = strlen(fname);
    if (len > 6 && !strncmp(fname + len - 6, ".gretl", 6)) {
	dirname[strlen(dirname) - 6] = '\0';
    } else {
	strcat(fname, ".gretl");
    }

#if SAVE_DEBUG
    fprintf(stderr, " save dirname = '%s'\n", dirname);
    fprintf(stderr, " current session.dirname = '%s'\n", session.dirname);
    fprintf(stderr, " doing chdir to '%s'\n", gretl_dotdir());
#endif

    /* note: paths below are relative to this */
    err = gretl_chdir(gretl_dotdir());
    if (err) {
	fprintf(stderr, " chdir to dotdir failed\n");
	gui_errmsg(err);
	return 1;
    }

    if (strcmp(dirname, session.dirname)) {
	/* have to rename the session directory */
	maybe_suspend_session_log();
	log_code = LOG_SAVE_AS;
	dirbak = gretl_strdup(session.dirname);
	if (gretl_isdir(dirname)) {
	    /* don't trip over a non-empty dir */
	    gretl_deltree(dirname);
	    gretl_error_clear();
	}
#ifdef G_OS_WIN32
	err = win32_rename_dir(session.dirname, dirname);
#else
	err = gretl_rename(session.dirname, dirname);
#endif
	if (err) {
	    fprintf(stderr, " failed to rename session dir: '%s' -> '%s'\n",
		    session.dirname, dirname);
	    gui_errmsg(err);
	} else {
	    fprintf(stderr, " renamed session dir OK\n");
	    strcpy(session.dirname, dirname);
	}
    }

    if (!err) {
	*datname = '\0';
	if (data_status) {
	    make_session_dataname(datname);
	} else {
	    strcpy(datname, "none");
	}
	err = write_session_xml(datname);
#if SAVE_DEBUG
	fprintf(stderr, " write_session_xml: err = %d\n", err);
#endif
    }

    if (!err && data_status) {
	err = real_save_session_dataset(datname);
#if SAVE_DEBUG
	fprintf(stderr, " real_save_session_dataset: err = %d\n", err);
#endif
	if (err) {
	    gui_errmsg(err);
	}
    }

    if (!err) {
	session_switch_log_location(log_code);
    }

    if (!err) {
	/* make zipfile containing session files */
	err = gretl_make_zipfile(fname, dirname);
	if (err) {
	    fprintf(stderr, " gretl_make_zipfile: err = %d\n", err);
	    gui_errmsg(err);
	} else {
	    mkfilelist(FILE_LIST_SESSION, fname, 0);
	    if (fname != sessionfile) {
		session_name_from_session_file(session.name, fname);
		strcpy(sessionfile, fname);
		data_status |= SESSION_DATA; /* FIXME? */
		set_sample_label(dataset);
	    }
	    mark_session_saved();
	}
    }

    if (dirbak != NULL) {
	if (err) {
	    /* restore original name on error */
	    strcpy(session.dirname, dirbak);
	}
	free(dirbak);
    }

    fprintf(stderr, "save_session; returning %d\n", err);

    return err;
}

void session_notes_callback (GtkWidget *w, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    const char *opts[] = {
	N_("Display notes on opening session file")
    };
    int active[] = {0};
    int resp;

    active[0] = session.show_notes;

    resp = checks_only_dialog("gretl", NULL, opts, 1,
			      active, 0, vwin->main);

    if (resp >= 0 && session.show_notes != active[0]) {
	session.show_notes = active[0];
	mark_session_changed();
    }
}

void save_session_callback (GtkAction *action)
{
    int as_is = 0;

    if (session_file_is_open() && action != NULL) {
	const gchar *s = gtk_action_get_name(action);

	if (!strcmp(s, "SaveSession")) {
	    as_is = 1;
	}
    }

    if (as_is) {
	save_session(sessionfile);
    } else {
	file_selector(SAVE_SESSION, FSEL_DATA_NONE, NULL);
    }
}

int save_session_commands (char *fname)
{
    FILE *fp = gretl_fopen(fname, "w");
    int err = 0;

    if (fp == NULL) {
	file_write_errbox(fname);
	err = E_FOPEN;
    } else {
	gchar *s = get_logfile_content(&err);

	if (err) {
	    gui_errmsg(err);
	} else {
	    fputs(s, fp);
	    g_free(s);
	}

	fclose(fp);
    }

    return err;
}

static char *model_cmd_str (MODEL *pmod)
{
    char *str = NULL;

    if (pmod->ci == MLE || pmod->ci == GMM || pmod->ncoeff > 10 ||
	pmod->list[0] > 10) {
	return NULL;
    }

    str = malloc(MAXLEN);
    if (str == NULL) {
	return NULL;
    }

    sprintf(str, "%s ", gretl_command_word(pmod->ci));

    if (pmod->ci == AR) {
        model_list_to_string(pmod->arinfo->arlist, str);
        strcat(str, "; ");
    }

    model_list_to_string(pmod->list, str);

    return str;
}

static gchar *graph_str (SESSION_GRAPH *graph)
{
    char tmp[MAXLEN];
    FILE *fp;
    gchar *buf = NULL;

    session_file_make_path(tmp, graph->fname, NULL);
    fp = gretl_fopen(tmp, "r");

    /* FIXME boxplots */

    if (fp != NULL) {
	char line[128], title[64];
	char xlabel[24], ylabel[24];
	int gottitle = 0;
	int gotxy = 0;

	while (fgets(line, sizeof line, fp)) {
	    if (strstr(line, "# timeseries") || strstr(line, "# frequency") ||
		!strncmp(line, "plot", 4)) {
		break;
	    } else if (sscanf(line, "set title '%63[^']", title) == 1) {
		gottitle = 1;
		break;
	    } else if (sscanf(line, "set xlabel '%23[^']", xlabel) == 1) {
		gotxy++;
	    } else if (sscanf(line, "set ylabel '%23[^']", ylabel) == 1) {
		gotxy++;
	    }
	}

	if (gottitle) {
	    buf = g_strdup(title);
	} else if (gotxy == 2) {
	    buf = g_strdup_printf("%s %s %s", ylabel, _("versus"), xlabel);
	}

	if (buf != NULL && !g_utf8_validate(buf, -1, NULL)) {
	    /* let's give up! */
	    g_free(buf);
	    buf = NULL;
	}

	fclose(fp);
    }

    return buf;
}

static int maybe_raise_object_window (gpointer data)
{
    windata_t *vwin = get_viewer_for_data(data);

    if (vwin != NULL) {
	gretl_viewer_present(vwin);
	return 1;
    } else {
	return 0;
    }
}

static int display_session_model (SESSION_MODEL *sm)
{
    DATASET *dset = dataset;
    PRN *prn;

    if (maybe_raise_object_window(sm->ptr)) {
	return 0;
    }

    if (sm->type != GRETL_OBJ_SYS && bufopen(&prn)) {
	return 1;
    }

    if (sm->type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) sm->ptr;

	if (pmod->submask == NULL && complex_subsampled()) {
	    dset = fetch_full_dataset();
	}
	printmodel(pmod, dset, OPT_NONE, prn);
	view_model(prn, pmod, sm->name);
    } else if (sm->type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) sm->ptr;

	gretl_VAR_print(var, dset, OPT_NONE, prn);
	view_buffer(prn, 78, 450, sm->name, var->ci, var);
    } else if (sm->type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) sm->ptr;

	edit_dialog(ESTIMATE, sm->name, NULL, NULL,
		    do_saved_eqn_system, sys,
		    VARCLICK_NONE, iconview);
    }

    return 0;
}

/* callback used in objectsave.c */

void session_model_callback (void *ptr, int action)
{
    SESSION_MODEL *mod = get_session_model_by_data(ptr);

    if (mod == NULL) {
	return;
    }

    if (action == OBJ_ACTION_SHOW) {
	display_session_model(mod);
    } else if (action == OBJ_ACTION_FREE) {
	delete_icon_for_data(mod);
	real_delete_model_from_session(mod);
    }
}

static void open_matrix (gui_obj *obj)
{
    user_var *u = (user_var *) obj->data;
    const char *name = user_var_get_name(u);

    edit_or_view_matrix(name, iconview);
}

static int is_dbnomics_bundle (const gretl_bundle *b)
{
    const char *s = gretl_bundle_get_creator((gretl_bundle *) b);

    return s != NULL && !strcmp(s, "dbnomics");
}

static void open_bundle (gui_obj *obj)
{
    user_var *u = (user_var *) obj->data;
    const char *name = user_var_get_name(u);
    gretl_bundle *b = user_var_get_value(u);
    int role = VIEW_BUNDLE;
    int done = 0;
    PRN *prn = NULL;

    if (maybe_raise_object_window(b)) {
	return;
    } else if (!gretl_bundle_has_content(b)) {
	warnbox(_("Bundle is empty"));
	return;
    } else if (bufopen(&prn)) {
	return;
    }

    done = try_exec_bundle_print_function(b, prn);

    if (!done) {
	/* nothing fancy, just show content */
	gretl_bundle_print(b, prn);
    } else if (is_dbnomics_bundle(b)) {
	role = VIEW_DBNOMICS;
    }

    view_buffer(prn, 80, 400, name, role, b);
}

static void open_gui_text (gui_obj *obj)
{
    SESSION_TEXT *text = (SESSION_TEXT *) obj->data;
    PRN *prn;

    prn = gretl_print_new_with_buffer(g_strdup(text->buf));

    if (prn != NULL) {
	view_buffer(prn, 80, 400, obj->name, INFO, NULL);
    }
}

static int real_delete_model_from_session (SESSION_MODEL *model)
{
    int nm = session.nmodels - 1;

    if (nm == 0) {
	free_session_model(session.models[0]);
	free(session.models);
	session.models = NULL;
    } else {
	SESSION_MODEL **mods;
	int i, j = 0;

	for (i=0; i<session.nmodels; i++) {
	    if (session.models[i]->ptr == model->ptr) {
		free_session_model(session.models[i]);
		j = i;
		break;
	    }
	}

	for (i=j; i<nm; i++) {
	    session.models[i] = session.models[i+1];
	}

	mods = myrealloc(session.models, nm * sizeof *mods);
	if (mods != NULL) {
	    session.models = mods;
	}
    }

    /* FIXME should delete the model's file in the session
       directory? */

    session.nmodels = nm;
    mark_session_changed();

    return 0;
}

static int real_delete_text_from_session (SESSION_TEXT *junk)
{
    int nt = session.ntexts;

    if (nt == 1) {
	free_session_text(session.texts[0]);
	free(session.texts);
	session.texts = NULL;
    } else {
	SESSION_TEXT **pptext;
	int i, j;

	pptext = mymalloc((nt - 1) * sizeof *pptext);
	if (pptext == NULL) {
	    return 1;
	}
	j = 0;
	for (i=0; i<nt; i++) {
	    if (session.texts[i] != junk) {
		pptext[j++] = session.texts[i];
	    } else {
		free_session_text(session.texts[i]);
	    }
	}
	free(session.texts);
	session.texts = pptext;
    }

    session.ntexts = nt - 1;
    mark_session_changed();

    return 0;
}

static void remove_session_graph_file (SESSION_GRAPH *graph)
{
    char fname[MAXLEN];

    gretl_chdir(gretl_dotdir());
    session_file_make_path(fname, graph->fname, NULL);
    gretl_remove(fname);

    if (graph->has_datafile) {
	gchar *datfile = g_strdup_printf("%s.dat", fname);

	gretl_remove(datfile);
	g_free(datfile);
    }
}

static int real_delete_graph_from_session (SESSION_GRAPH *junk)
{
    int ng = session.ngraphs;

    if (in_graph_page(junk->fname)) {
	graph_page_delete_file(junk->fname);
    }

    if (ng == 1) {
	remove_session_graph_file(session.graphs[0]);
	free(session.graphs[0]);
	free(session.graphs);
	session.graphs = NULL;
    } else {
	SESSION_GRAPH **ppgr;
	int i, j, done = 0;

	for (i=0; i<ng && !done; i++) {
	    if (!strcmp(session.graphs[i]->name, junk->name)) {
		remove_session_graph_file(session.graphs[i]);
		free(session.graphs[i]);
		for (j=i; j<ng-1; j++) {
		    session.graphs[j] = session.graphs[j+1];
		}
		done = 1;
	    }
	}

	if (done) {
	    ppgr = myrealloc(session.graphs, (ng - 1) * sizeof *ppgr);
	    if (ppgr == NULL) {
		return 1;
	    }
	    session.graphs = ppgr;
	}
    }

    session.ngraphs = ng - 1;
    mark_session_changed();

    return 0;
}

static int delete_session_object (gui_obj *obj)
{
    if (obj->sort == GRETL_OBJ_EQN || obj->sort == GRETL_OBJ_VAR ||
	obj->sort == GRETL_OBJ_SYS) {
	real_delete_model_from_session(obj->data);
    } else if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	real_delete_graph_from_session(obj->data);
    } else if (obj->sort == GRETL_OBJ_TEXT) {
	real_delete_text_from_session(obj->data);
    } else if (obj->sort == GRETL_OBJ_MATRIX ||
	       obj->sort == GRETL_OBJ_BUNDLE) {
	user_var_delete(obj->data);
    }

    session_delete_icon(obj);
    mark_session_changed();

    return 0;
}

/* run a sanity check on deleting a session object (e.g.
   model, graph) before proceeding */

static void maybe_delete_session_object (gui_obj *obj)
{
    GtkWidget *busywin = NULL;

    if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	SESSION_GRAPH *graph = (SESSION_GRAPH *) obj->data;

	busywin = get_window_for_plot(graph);
	if (busywin == NULL) {
	    char fullname[MAXLEN];

	    session_file_make_path(fullname, graph->fname, NULL);
	    busywin = vwin_toplevel(get_editor_for_file(fullname));
	    if (busywin == NULL) {
		busywin = get_viewer_for_plot(fullname);
	    }
	}
    } else {
	gpointer p = NULL;

	if (obj->sort == GRETL_OBJ_EQN || obj->sort == GRETL_OBJ_SYS ||
	    obj->sort == GRETL_OBJ_VAR) {
	    SESSION_MODEL *mod = obj->data;

	    p = mod->ptr;
	} else if (obj->sort == GRETL_OBJ_BUNDLE) {
	    p = user_var_get_value((user_var *) obj->data);
	} else {
	    p = obj->data;
	}

	if (p != NULL) {
	    if (obj->sort == GRETL_OBJ_MATRIX) {
		busywin = get_window_for_data(p);
	    } else {
		busywin = vwin_toplevel(get_viewer_for_data(p));
	    }
	}
    }

    if (busywin != NULL) {
	gtk_window_present(GTK_WINDOW(busywin));
	warnbox_printf(_("%s: please close this object's window first"),
		       obj->name);
    } else {
	gchar *msg = NULL;

	if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	    SESSION_GRAPH *graph = (SESSION_GRAPH *) obj->data;

	    if (in_graph_page(graph->fname)) {
		msg = g_strdup_printf(_("Really delete %s?\n(This will remove it from "
					"the Graph page.)"), obj->name);
	    }
	}
	if (msg == NULL) {
	    msg = g_strdup_printf(_("Really delete %s?"), obj->name);
	}
	if (yes_no_dialog(_("gretl: delete"), msg, iconview) == GRETL_YES) {
	    delete_session_object(obj);
	}
	g_free(msg);
    }
}

static gui_obj *get_gui_obj_by_data (void *targ)
{
    GList *mylist = g_list_first(iconlist);
    gui_obj *obj = NULL;

    while (mylist != NULL) {
	obj = (gui_obj *) mylist->data;
	if (obj->data == targ) {
	    return obj;
	}
	mylist = mylist->next;
    }

    return NULL;
}

static GretlObjType get_obj_type (GretlType type)
{
    if (type == GRETL_TYPE_MATRIX) {
	return GRETL_OBJ_MATRIX;
    } else if (type == GRETL_TYPE_BUNDLE) {
	return GRETL_OBJ_BUNDLE;
    } else {
	return GRETL_OBJ_NULL;
    }
}

/* called from DELEET case in library.c */

int session_user_var_destroy_by_name (const char *name,
				      GretlObjType type)
{
    user_var *u = get_user_var_by_name(name);
    int err;

    if (u == NULL) {
	err = E_UNKVAR;
    } else {
	maybe_close_window_for_user_var(u, type);
	if (iconlist != NULL) {
	    gui_obj *obj = get_gui_obj_by_data(u);

	    session_delete_icon(obj);
	}
	err = user_var_delete(u);
    }

    return err;
}

static void rename_session_graph (SESSION_GRAPH *graph, const char *newname)
{
    int i;

    for (i=0; i<session.ngraphs; i++) {
	if (!strcmp(session.graphs[i]->name, graph->name)) {
	    session.graphs[i]->name[0] = '\0';
	    strncat(session.graphs[i]->name, newname, MAXSAVENAME - 1);
	    break;
	}
    }
}

static void maybe_sync_model_window_name (SESSION_MODEL *sm)
{
    windata_t *vwin = get_viewer_for_data(sm->ptr);

    if (vwin != NULL) {
	gchar *title = g_strdup_printf("gretl: %s", sm->name);

	gretl_viewer_set_title(vwin, title);
	g_free(title);
    }
}

static int rename_session_object (gui_obj *obj, const char *newname,
				  GtkWidget *parent)
{
    int err_shown = 0;
    int err = 0;

    if (obj->sort == GRETL_OBJ_EQN || obj->sort == GRETL_OBJ_SYS ||
	obj->sort == GRETL_OBJ_VAR) {
	SESSION_MODEL *sm;

	sm = get_session_model_by_name(newname);
	if (sm != NULL) {
	    err = 1;
	} else {
	    sm = obj->data;
	    gretl_object_rename(sm->ptr, sm->type, newname);
	    *sm->name = '\0';
	    strncat(sm->name, newname, MAXSAVENAME - 1);
	    maybe_sync_model_window_name(sm);
	}
    } else if (obj->sort == GRETL_OBJ_GRAPH ||
	       obj->sort == GRETL_OBJ_PLOT) {
	SESSION_GRAPH *sg;

	sg = get_session_graph_by_name(newname);
	if (sg != NULL) {
	    err = 1;
	} else {
	    sg = obj->data;
	    rename_session_graph(sg, newname);
	}
    } else if (obj->sort == GRETL_OBJ_MATRIX ||
	       obj->sort == GRETL_OBJ_BUNDLE) {
	GretlType type = (obj->sort == GRETL_OBJ_MATRIX)?
	    GRETL_TYPE_MATRIX : GRETL_TYPE_BUNDLE;

	err = gui_validate_varname_strict(newname, type, parent);
	if (err) {
	    err_shown = 1;
	} else {
	    user_var_set_name(obj->data, newname);
	}
    } else if (obj->sort == GRETL_OBJ_TEXT) {
	SESSION_TEXT *st;

	st = get_session_text_by_name(newname);
	if (st != NULL) {
	    err = 1;
	} else {
	    st = obj->data;
	    *st->name = '\0';
	    strncat(st->name, newname, MAXSAVENAME - 1);
	}
    }

    if (err && !err_shown) {
	errbox_printf(_("'%s': there is already an object of this name"),
		      newname);
    } else if (!err) {
	free(obj->name);
	obj->name = g_strdup(newname);
    }

    return err;
}

static int copy_session_object (gui_obj *obj, const char *cpyname)
{
    void *oldp = NULL;
    void *p = NULL;
    int ptype = 0;
    int err = 0;

    /* Only graphs and matrices are supported for copying */

    if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	oldp = get_session_graph_by_name(cpyname);
    } else if (obj->sort == GRETL_OBJ_MATRIX) {
	oldp = get_user_var_by_name(cpyname);
    }

    if (oldp != NULL) {
	errbox_printf(_("'%s': there is already an object of this name"),
		      cpyname);
    } else if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	SESSION_GRAPH *g0 = obj->data;

	errno = 0;
	if (!session_dir_ok()) {
	    err = 1;
	} else {
	    char fname1[MAXSAVENAME];
	    char path0[FILENAME_MAX];
	    char path1[FILENAME_MAX];

	    err = gretl_chdir(gretl_dotdir());
	    if (!err) {
		make_graph_filename(fname1);
		session_file_make_path(path0, g0->fname, NULL);
		session_file_make_path(path1, fname1, NULL);
		err = copyfile(path0, path1);
	    }
	    if (!err) {
		p = session_append_graph(cpyname, fname1, g0->type);
		ptype = g0->type;
		err = (p == NULL);
	    }
	}
    } else if (obj->sort == GRETL_OBJ_MATRIX) {
	user_var *u = obj->data;

	err = copy_matrix_as(user_var_get_value(u), cpyname, 0);
	if (!err) {
	    p = get_user_var_by_name(cpyname);
	    ptype = GRETL_OBJ_MATRIX;
	    err = (p == NULL);
	}
    }

    if (!err) {
	session_add_icon(p, ptype, ICON_ADD_SINGLE);
    }

    return err;
}

static void copy_object_callback (GtkWidget *widget, dialog_t *dlg)
{
    gui_obj *obj = (gui_obj *) edit_dialog_get_data(dlg);
    const gchar *cpyname;
    int err = 0;

    cpyname = edit_dialog_get_text(dlg);

    if (cpyname != NULL && *cpyname != '\0') {
	err = copy_session_object(obj, cpyname);
	if (err) {
	    errbox(_("Failed to copy object"));
	} else {
	    mark_session_changed();
	}
    }

    if (!err) {
	edit_dialog_close(dlg);
    }
}

static void copy_object_dialog (gui_obj *obj)
{
    int maxlen = MAXSAVENAME - 1;
    gchar *tmp;

    if (obj->sort != GRETL_OBJ_MATRIX &&
	obj->sort != GRETL_OBJ_GRAPH &&
	obj->sort != GRETL_OBJ_PLOT) {
	dummy_call();
    }

    if (obj->sort == GRETL_OBJ_MATRIX || obj->sort == GRETL_OBJ_BUNDLE) {
	maxlen = VNAMELEN - 1;
    }

    tmp = g_strdup_printf(_("Enter new name\n(max. %d characters)"),
			  maxlen);
    edit_dialog(0, _("gretl: copy object"), tmp, obj->name,
		copy_object_callback, obj,
		VARCLICK_NONE, iconview);
    g_free(tmp);
}

static void rename_object_callback (GtkWidget *widget, dialog_t *dlg)
{
    gui_obj *obj = (gui_obj *) edit_dialog_get_data(dlg);
    const gchar *newname;
    int err = 0;

    newname = edit_dialog_get_text(dlg);

    if (newname != NULL && *newname != '\0' &&
	strcmp(newname, obj->name)) {
	GtkWidget *parent = edit_dialog_get_window(dlg);
	gchar str[2*SHOWNAMELEN];

	err = rename_session_object(obj, newname, parent);
	if (!err) {
	    make_short_label_string(str, obj->name);
	    gtk_label_set_text(GTK_LABEL(obj->label), str);
	    mark_session_changed();
	}
    }

    if (!err) {
	edit_dialog_close(dlg);
    }
}

static void rename_object_dialog (gui_obj *obj)
{
    int maxlen = MAXSAVENAME - 1;
    gchar *tmp;

    if (obj->sort == GRETL_OBJ_MATRIX || obj->sort == GRETL_OBJ_BUNDLE) {
	maxlen = VNAMELEN - 1;
    }

    tmp = g_strdup_printf(_("Enter new name\n(max. %d characters)"),
			  maxlen);
    edit_dialog(0, _("gretl: rename object"), tmp, obj->name,
		rename_object_callback, obj,
		VARCLICK_NONE, iconview);
    g_free(tmp);
}

void delete_text_from_session (void *p)
{
    SESSION_TEXT *text = (SESSION_TEXT *) p;
    gui_obj *obj;

    if (text == NULL) return;

    real_delete_text_from_session(text);

    obj = get_gui_obj_by_data((void *) text);
    if (obj != NULL) {
	session_delete_icon(obj);
    }
}

void display_saved_text (void *p)
{
    SESSION_TEXT *text = (SESSION_TEXT *) p;
    PRN *prn;

    if (text == NULL) return;

    prn = gretl_print_new_with_buffer(g_strdup(text->buf));
    if (prn != NULL) {
	view_buffer(prn, 80, 400, text->name, INFO, NULL);
    }
}

static void session_view_init (void)
{
    iconlist = NULL;
    icon_table = NULL;
    iconview_width = 0;
    iconview_cols = ICONVIEW_MIN_COLS;
}

static void gui_obj_destroy (gui_obj *obj, gpointer p)
{
    if (obj != NULL) {
#if SESSION_DEBUG
	fprintf(stderr, "freeing obj at %p (%s)\n", (void *) obj, obj->name);
#endif
	if (obj->name != NULL) {
	    g_free(obj->name);
	}
	g_object_unref(obj->icon);
	g_object_unref(obj->label);
	free(obj);
    }
}

static void session_view_free (GtkWidget *w, gpointer data)
{
    iconview = NULL;

    g_list_foreach(iconlist, (GFunc) gui_obj_destroy, NULL);

    g_list_free(iconlist);
    iconlist = NULL;
}

static void session_delete_icon (gui_obj *obj)
{
    if (obj == NULL) return;

    if (obj->icon != NULL && GTK_IS_WIDGET(obj->icon)) {
	gtk_container_remove(GTK_CONTAINER(icon_table), obj->icon);
    }
    if (obj->label != NULL && GTK_IS_WIDGET(obj->label)) {
	gtk_container_remove(GTK_CONTAINER(icon_table), obj->label);
    }

    iconlist = g_list_first(iconlist);
    iconlist = g_list_remove(iconlist, obj);

    gui_obj_destroy(obj, NULL);

    if (iconlist == NULL) {
	fprintf(stderr, "Bad: iconlist has gone NULL\n");
    }
}

/* apparatus for getting a white background */

#if GTK_MAJOR_VERSION >= 3

static void white_bg_style (GtkWidget *widget, gpointer data)
{
    static GdkRGBA rgbw = {1, 1, 1, 1};
    static GdkRGBA rgbb;
    static int done;

    gtk_widget_override_background_color(widget,
					 GTK_STATE_FLAG_NORMAL,
					 &rgbw);
    if (!done) {
	gdk_rgba_parse(&rgbb, "#4a90d9");
    }
    gtk_widget_override_background_color(widget,
					 GTK_STATE_FLAG_PRELIGHT,
					 &rgbb);
}

#else

static GdkColor *get_white (void)
{
    GdkColormap *cmap;
    GdkColor *white;

    white = mymalloc(sizeof *white);
    if (white == NULL) return NULL;

    cmap = gdk_colormap_get_system();
    gdk_color_parse("white", white);
    gdk_colormap_alloc_color(cmap, white, FALSE, TRUE);

    return white;
}

static void white_bg_style (GtkWidget *widget, gpointer data)
{
    static GdkColor *white;

    if (white == NULL) {
	white = get_white();
    }

    gtk_widget_modify_bg(widget, GTK_STATE_NORMAL, white);
}

#endif

static void real_pack_icon (gui_obj *obj, int row, int col)
{
    obj->row = row;
    obj->col = col;

    gtk_table_attach(GTK_TABLE(icon_table), obj->icon,
		     col, col + 1, row, row + 1,
		     GTK_EXPAND, GTK_FILL, 5, 5);

    gtk_widget_show(obj->icon);
    white_bg_style(obj->icon, NULL);

    gtk_table_attach(GTK_TABLE(icon_table), obj->label,
		     col, col + 1, row + 1, row + 2,
		     GTK_EXPAND, GTK_FILL, 5, 5);

    gtk_widget_show(obj->label);
}

static void pack_single_icon (gui_obj *obj)
{
    int row, col;
    gui_obj *last;

    iconlist = g_list_last(iconlist);
    last = iconlist->data;
    row = last->row;
    col = last->col;

    iconlist = g_list_append(iconlist, obj);

    col++;
    if (col > 0 && (col % iconview_cols == 0)) {
	col = 0;
	row += 2;
	gtk_table_resize(GTK_TABLE(icon_table), 2 * row, iconview_cols);
    }

    real_pack_icon(obj, row, col);
}

/* returns the number of rows used, counting 1 per icon + label */

static int batch_pack_icons (void)
{
    GList *mylist = g_list_first(iconlist);
    gui_obj *obj;
    int row = 0, col = 0;

    while (mylist != NULL) {
	obj = (gui_obj *) mylist->data;
	real_pack_icon(obj, row, col);
	col++;
	if (col > 0 && (col % iconview_cols == 0)) {
	    col = 0;
	    row += 2;
	    gtk_table_resize(GTK_TABLE(icon_table), 2 * row, iconview_cols);
	}
	if (mylist->next == NULL) {
	    break;
	} else {
	    mylist = mylist->next;
	}
    }

    return row / 2;
}

static int gui_user_var_callback (const char *name, GretlType type,
				  int action)
{
    GretlObjType otype = get_obj_type(type);
    user_var *u;
    int err = 0;

    if (action == UVAR_ADD) {
	/* variable has been added, GUI sync wanted */
	if (iconview != NULL) {
	    u = get_user_var_of_type_by_name(name, type);
	    if (u != NULL) {
		session_add_icon(u, otype, ICON_ADD_SINGLE);
	    }
	} else if (type != GRETL_TYPE_MATRIX && autoicon_on()) {
	    /* auto-open icon view unless we added a matrix */
	    auto_view_session();
	}
	if (iconview != NULL && waiting_for_output()) {
	    gtk_widget_set_sensitive(iconview, FALSE);
	}
	mark_session_changed();
    } else if (action == UVAR_DELETE) {
	/* variable not yet deleted (deferred to GUI) */
	u = get_user_var_of_type_by_name(name, type);
	if (u == NULL) {
	    err = E_UNKVAR;
	} else {
	    maybe_close_window_for_user_var(u, otype);
	    if (iconlist != NULL) {
		session_delete_icon(get_gui_obj_by_data(u));
	    }
	    err = user_var_delete(u);
	}
    }

    return err;
}

void maybe_sensitize_iconview (void)
{
    if (iconview != NULL && !gtk_widget_is_sensitive(iconview)) {
	gtk_widget_set_sensitive(iconview, TRUE);
    }
}

int have_session_objects (void)
{
    int n = data_status;

    /* Note: the following lines enable the icon view even when
       there's no dataset present.
    */

    if (n == 0) {
	n = session.ngraphs + session.ntexts;
    }

    if (n == 0) {
	n = n_user_matrices() + n_user_bundles() + n_user_scalars();
    }

    return n > 0;
}

static void add_user_var_icon (gpointer data, gpointer intp)
{
    session_add_icon(data, GPOINTER_TO_INT(intp), ICON_ADD_BATCH);
}

/* returns the number of rows of (icon + label) */

static int add_all_icons (void)
{
    int show_graph_page = check_for_program(latex);
    GList *list = NULL;
    int i;

    active_object = NULL;

    if (data_status) {
	session_add_icon(NULL, GRETL_OBJ_INFO,    ICON_ADD_BATCH);  /* data info */
	session_add_icon(NULL, GRETL_OBJ_DSET,    ICON_ADD_BATCH);  /* data file */
	session_add_icon(NULL, GRETL_OBJ_STATS,   ICON_ADD_BATCH);  /* summary stats */
	session_add_icon(NULL, GRETL_OBJ_CORR,    ICON_ADD_BATCH);  /* correlation matrix */
	session_add_icon(NULL, GRETL_OBJ_MODTAB,  ICON_ADD_BATCH);  /* model table */
    }

    /* standard icons that don't really require a dataset in place */
    session_add_icon(NULL, GRETL_OBJ_SCALARS, ICON_ADD_BATCH);  /* scalars */
    session_add_icon(NULL, GRETL_OBJ_NOTES,   ICON_ADD_BATCH);  /* session notes */
    if (show_graph_page) {
	session_add_icon(NULL, GRETL_OBJ_GPAGE, ICON_ADD_BATCH); /* graph page */
    }

    for (i=0; i<session.nmodels; i++) {
#if SESSION_DEBUG
	fprintf(stderr, "adding session.models[%d] (type %d) to view\n", i,
		session.models[i]->type);
#endif
	session_add_icon(session.models[i], session.models[i]->type,
			 ICON_ADD_BATCH);
    }

    for (i=0; i<session.ngraphs; i++) {
#if SESSION_DEBUG
	fprintf(stderr, "adding session.graphs[%d] (type %d) to view\n", i,
		session.graphs[i]->type);
#endif
	session_add_icon(session.graphs[i], session.graphs[i]->type,
			 ICON_ADD_BATCH);
    }

    for (i=0; i<session.ntexts; i++) {
#if SESSION_DEBUG
	fprintf(stderr, "adding session.texts[%d] to view\n", i);
#endif
	session_add_icon(session.texts[i], GRETL_OBJ_TEXT, ICON_ADD_BATCH);
    }

    list = user_var_list_for_type(GRETL_TYPE_MATRIX);
    g_list_foreach(list, add_user_var_icon, GINT_TO_POINTER(GRETL_OBJ_MATRIX));
    g_list_free(list);

    list = user_var_list_for_type(GRETL_TYPE_BUNDLE);
    g_list_foreach(list, add_user_var_icon, GINT_TO_POINTER(GRETL_OBJ_BUNDLE));
    g_list_free(list);

    return batch_pack_icons();
}

static void undisplay_icon (gui_obj *obj, gpointer p)
{
    if (obj->icon && GTK_IS_WIDGET(obj->icon))
	gtk_container_remove(GTK_CONTAINER(icon_table), obj->icon);
    if (obj->label && GTK_IS_WIDGET(obj->label))
	gtk_container_remove(GTK_CONTAINER(icon_table), obj->label);
}

static void rearrange_icons (void)
{
    iconlist = g_list_first(iconlist);
    g_list_foreach(iconlist, (GFunc) undisplay_icon, NULL);
    batch_pack_icons();
}

static gint catch_iconview_key (GtkWidget *w, GdkEventKey *key,
				gpointer p)
{
#ifdef __APPLE__
    if (key->keyval == GDK_w && cmd_key(key)) {
	gtk_widget_destroy(w);
    }
#endif
    /* 'q' quits iconview */
    if (key->keyval == GDK_q) {
        gtk_widget_destroy(w);
    }

    return FALSE;
}

static void mtab_item_set_sensitivity (gui_obj *obj)
{
    SESSION_MODEL *model = obj->data;
    int added = in_model_table(model->ptr);

    if (mtab_add_item != NULL) {
	gtk_widget_set_sensitive(mtab_add_item, !added);
    }
}

static void gpage_item_set_sensitivity (gui_obj *obj)
{
    SESSION_GRAPH *graph = obj->data;
    int added = in_graph_page(graph->fname);

    if (gpage_add_item != NULL) {
	gtk_widget_set_sensitive(gpage_add_item, !added);
    }
}

static void copy_item_set_sensitivity (gui_obj *obj)
{
    SESSION_GRAPH *graph = obj->data;

    if (graph_copy_item != NULL) {
	gtk_widget_set_sensitive(graph_copy_item, !graph->has_datafile);
    }
}


static void object_popup_show (gui_obj *obj, GdkEventButton *event)
{
    GtkWidget *w = NULL;

    active_object = obj;

    switch (obj->sort) {
    case GRETL_OBJ_EQN:
	w = model_popup;
	mtab_item_set_sensitivity(obj);
	break;
    case GRETL_OBJ_MODTAB:
	w = model_table_popup;
	break;
    case GRETL_OBJ_GPAGE:
	w = graph_page_popup;
	break;
    case GRETL_OBJ_VAR:
    case GRETL_OBJ_SYS:
    case GRETL_OBJ_TEXT:
	w = generic_popup;
	break;
    case GRETL_OBJ_GRAPH:
    case GRETL_OBJ_PLOT:
	gpage_item_set_sensitivity(obj);
	copy_item_set_sensitivity(obj);
	w = graph_popup;
	break;
    case GRETL_OBJ_DSET:
	w = data_popup;
	break;
    case GRETL_OBJ_SCALARS:
	w = scalars_popup;
	break;
    case GRETL_OBJ_INFO:
	w = info_popup;
	break;
    case GRETL_OBJ_MATRIX:
	w = matrix_popup;
	break;
    case GRETL_OBJ_BUNDLE:
	w = bundle_popup;
	break;
    default:
	break;
    }

    gtk_menu_popup(GTK_MENU(w), NULL, NULL, NULL, NULL,
		   event->button, event->time);
}

static void display_model_table_wrapper (void)
{
    display_model_table(1);
}

static void graph_page_save_wrapper (void)
{
    if (graph_page_get_n_graphs() == 0) {
	warnbox(_("The graph page is empty"));
    } else {
	file_selector(SAVE_TEX, FSEL_DATA_NONE, NULL);
    }
}

static gboolean session_icon_click (GtkWidget *icon,
				    GdkEventButton *event,
				    gpointer data)
{
    gui_obj *obj = (gui_obj *) data;

    if (event->type == GDK_2BUTTON_PRESS) {
	switch (obj->sort) {
	case GRETL_OBJ_EQN:
	case GRETL_OBJ_VAR:
	case GRETL_OBJ_SYS:
	    display_session_model(obj->data);
	    break;
	case GRETL_OBJ_GRAPH:
	case GRETL_OBJ_PLOT:
	    open_gui_graph(obj);
	    break;
	case GRETL_OBJ_TEXT:
	    open_gui_text(obj);
	    break;
	case GRETL_OBJ_DSET:
	    show_spreadsheet(SHEET_EDIT_DATASET);
	    break;
	case GRETL_OBJ_SCALARS:
	    edit_scalars();
	    break;
	case GRETL_OBJ_MATRIX:
	    open_matrix(obj);
	    break;
	case GRETL_OBJ_BUNDLE:
	    open_bundle(obj);
	    break;
	case GRETL_OBJ_INFO:
	    dataset_info();
	    break;
	case GRETL_OBJ_NOTES:
	    edit_session_notes();
	    break;
	case GRETL_OBJ_MODTAB:
	    display_model_table_wrapper();
	    break;
	case GRETL_OBJ_GPAGE:
	    display_graph_page(iconview);
	    break;
	case GRETL_OBJ_CORR:
	    do_menu_op(ALL_CORR, NULL, OPT_NONE, iconview);
	    break;
	case GRETL_OBJ_STATS:
	    do_menu_op(ALL_SUMMARY, NULL, OPT_NONE, iconview);
	    break;
	}
	return TRUE;
    } else {
	if (right_click(event)) {
	    if (obj->sort == GRETL_OBJ_EQN  || obj->sort == GRETL_OBJ_GRAPH ||
		obj->sort == GRETL_OBJ_TEXT || obj->sort == GRETL_OBJ_DSET ||
		obj->sort == GRETL_OBJ_INFO || obj->sort == GRETL_OBJ_GPAGE ||
		obj->sort == GRETL_OBJ_PLOT || obj->sort == GRETL_OBJ_MODTAB ||
		obj->sort == GRETL_OBJ_VAR  || obj->sort == GRETL_OBJ_SYS ||
		obj->sort == GRETL_OBJ_MATRIX || obj->sort == GRETL_OBJ_BUNDLE ||
		obj->sort == GRETL_OBJ_SCALARS) {
		object_popup_show(obj, event);
	    }
	    return TRUE;
	}
    }

    return FALSE;
}

static gboolean session_view_click (GtkWidget *widget,
				    GdkEventButton *event,
				    gpointer data)
{
    if (!in_icon) {
	if (right_click(event)) {
	    /* right-click on iconview background */
	    gtk_menu_popup(GTK_MENU(global_popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    return TRUE;
	}
    }

    return FALSE;
}

static void global_popup_callback (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (!strcmp(item, _("Save session"))) {
	if (sessionfile[0]) {
	    if (session.status & SESSION_CHANGED) {
		save_session(NULL);
	    } else {
		mark_session_saved();
	    }
	} else {
	    file_selector(SAVE_SESSION, FSEL_DATA_NONE, NULL);
	}
    } else if (!strcmp(item, _("Arrange icons"))) {
	rearrange_icons();
    } else if (!strcmp(item, _("Add matrix..."))) {
	gui_new_matrix(iconview);
    } else if (!strcmp(item, _("Windows"))) {
	window_list_popup(widget, NULL, iconview);
    } else if (!strcmp(item, _("Close window"))) {
	gtk_widget_destroy(iconview);
    }
}

static void info_popup_callback (GtkWidget *widget, gpointer data)
{
    dataset_info();
}

static void matrix_popup_callback (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;
    gui_obj *obj = active_object;
    user_var *u = (user_var *) obj->data;
    const char *name = user_var_get_name(u);
    gretl_matrix *m;

    if (!strcmp(item, _("View"))) {
	PRN *prn;

	m = user_var_get_value(u);
	if (m != NULL && bufopen(&prn) == 0) {
	    gretl_matrix_print_to_prn(m, name, prn);
	    view_buffer(prn, 78, 400, name, PRINT, NULL);
	}
    } else if (!strcmp(item, _("Edit"))) {
	edit_or_view_matrix(name, iconview);
    } else if (!strcmp(item, _("Properties"))) {
	m = user_var_get_value(u);
	view_matrix_properties(m, name);
    } else if (!strcmp(item, _("Copy as CSV..."))) {
	m = user_var_get_value(u);
	if (gretl_is_null_matrix(m)) {
	    warnbox("matrix is null");
	} else {
	    matrix_to_clipboard_as_csv(m, iconview);
	}
    } else if (!strcmp(item, _("Rename"))) {
	rename_object_dialog(obj);
    } else if (!strcmp(item, _("Delete"))) {
	maybe_delete_session_object(obj);
    } else if (!strcmp(item, _("Copy"))) {
	copy_object_dialog(obj);
    }
}

static void bundle_popup_callback (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;
    gui_obj *obj = active_object;

    if (!strcmp(item, _("View"))) {
	open_bundle(obj);
    } else if (!strcmp(item, _("Rename"))) {
	rename_object_dialog(obj);
    } else if (!strcmp(item, _("Delete"))) {
	maybe_delete_session_object(obj);
    }
}

static void data_popup_callback (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (!strcmp(item, _("Edit"))) {
	show_spreadsheet(SHEET_EDIT_DATASET);
    } else if (!strcmp(item, _("Export as CSV..."))) {
	file_save(mdata, EXPORT_CSV);
    } else if (!strcmp(item, _("Copy as CSV..."))) {
	csv_to_clipboard(iconview);
    }
}

static void scalars_popup_callback (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (!strcmp(item, _("Edit"))) {
	edit_scalars();
    } else if (!strcmp(item, _("Copy as CSV..."))) {
	scalars_to_clipboard_as_csv(iconview);
    }
}

static gchar *object_get_window_title (gui_obj *obj)
{
    gchar *title = NULL;

    if (obj != NULL) {
	title = g_strdup_printf("gretl: %s", obj->name);
    }

    return title;
}

static void object_popup_callback (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;
    gui_obj *obj = active_object;

    if (!strcmp(item, _("Display"))) {
	if (obj->sort == GRETL_OBJ_EQN ||
	    obj->sort == GRETL_OBJ_VAR ||
	    obj->sort == GRETL_OBJ_SYS) {
	    display_session_model(obj->data);
	} else if (obj->sort == GRETL_OBJ_TEXT) {
	    open_gui_text(obj);
	} else if (obj->sort == GRETL_OBJ_MODTAB) {
	    display_model_table_wrapper();
	} else if (obj->sort == GRETL_OBJ_GPAGE) {
	    display_graph_page(iconview);
	} else if (obj->sort == GRETL_OBJ_GRAPH ||
		   obj->sort == GRETL_OBJ_PLOT) {
	    open_gui_graph(obj);
	}
    } else if (!strcmp(item, _("Edit plot commands"))) {
	if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	    SESSION_GRAPH *graph = (SESSION_GRAPH *) obj->data;
	    char fullname[MAXLEN];
	    gchar *title;
	    windata_t *vwin;
	    int err;

	    err = gretl_chdir(gretl_dotdir());
	    if (err) {
		gui_errmsg(err);
	    } else {
		session_file_make_path(fullname, graph->fname, NULL);
		/* the following handles error message if needed */
		err = remove_png_term_from_plot_by_name(fullname);
	    }
	    if (!err) {
		title = object_get_window_title(obj);
		vwin = view_file_with_title(fullname, 1, 0, 78, 400,
					    EDIT_GP, title);
		g_free(title);
		/* add flag so we can mark the session as modified
		   if the plot file is changed */
		vwin->flags |= VWIN_SESSION_GRAPH;
		if (graph->has_datafile) {
		    vwin->data = graph;
		}
	    }
	}
    } else if (!strcmp(item, _("Add to graph page"))) {
	if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	    SESSION_GRAPH *graph = (SESSION_GRAPH *) obj->data;

	    if (!in_graph_page(graph->fname)) {
		graph_page_add_file(graph->fname);
	    }
	}
    } else if (!strcmp(item, _("Rename"))) {
	rename_object_dialog(obj);
    } else if (!strcmp(item, _("Delete"))) {
	/* note: "Delete" = "Clear" in some translations */
	if (obj->sort == GRETL_OBJ_MODTAB) {
	    clear_model_table(0, NULL);
	} else if (obj->sort == GRETL_OBJ_GPAGE) {
	    clear_graph_page(0);
	} else {
	    maybe_delete_session_object(obj);
	}
    } else if (!strcmp(item, _("Add to model table"))) {
	if (obj->sort == GRETL_OBJ_EQN) {
	    SESSION_MODEL *mod = (SESSION_MODEL *) obj->data;

	    add_to_model_table(mod->ptr, MODEL_ADD_FROM_MENU, 0, NULL);
	}
    } else if (!strcmp(item, _("Clear"))) {
	if (obj->sort == GRETL_OBJ_MODTAB) {
	    clear_model_table(0, NULL);
	} else if (obj->sort == GRETL_OBJ_GPAGE) {
	    clear_graph_page(0);
	}
    } else if (!strcmp(item, _("Help"))) {
	if (obj->sort == GRETL_OBJ_MODTAB) {
	    show_gui_help(MODELTAB);
	} else if (obj->sort == GRETL_OBJ_GPAGE) {
	    show_gui_help(GRAPHPG);
	}
    } else if (!strcmp(item, _("Save as TeX..."))) {
	if (obj->sort == GRETL_OBJ_GPAGE) {
	    graph_page_save_wrapper();
	}
    } else if (!strcmp(item, _("Copy"))) {
	if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	    copy_object_dialog(obj);
	}
    }
}

static gboolean icon_entered (GtkWidget *icon, GdkEventCrossing *event,
			      gui_obj *obj)
{
#if GTK_MAJOR_VERSION == 3
    gtk_widget_set_state_flags(icon, GTK_STATE_FLAG_PRELIGHT, FALSE);
#else
    gtk_widget_set_state(icon, GTK_STATE_SELECTED);
#endif
    in_icon = 1;

    return FALSE;
}

static gboolean icon_left (GtkWidget *icon, GdkEventCrossing *event,
			   gui_obj *obj)
{
    gtk_widget_set_state(icon, GTK_STATE_NORMAL);
    in_icon = 0;

    return FALSE;
}

static MODEL *drag_model_src;

static void
session_data_received (GtkWidget *widget,
		       GdkDragContext *context,
		       gint x,
		       gint y,
		       GtkSelectionData *data,
		       guint info,
		       guint time,
		       gpointer p)
{
    const guchar *seldata = NULL;

    if (data != NULL) {
	seldata = gtk_selection_data_get_data(data);
    }

    if (info == GRETL_MODEL_PTR && seldata != NULL) {
	MODEL *pmod = drag_model_src;

	if (pmod != NULL) {
	    add_to_model_table(pmod, MODEL_ADD_BY_DRAG, 0, NULL);
	    drag_model_src = NULL;
	}
    } else if (info == GRETL_GRAPH_FILE && seldata != NULL) {
	gchar *fname = (gchar *) seldata;

	if (fname != NULL) {
	    graph_page_add_file(fname);
	}
    }
}

static void session_drag_setup (gui_obj *obj)
{
    GtkWidget *w = GTK_WIDGET(obj->icon);
    GtkTargetEntry *targ;

    if (obj->sort == GRETL_OBJ_MODTAB) {
	targ = &gretl_drag_targets[GRETL_MODEL_PTR];
    } else {
	targ = &gretl_drag_targets[GRETL_GRAPH_FILE];
    }

    gtk_drag_dest_set(w,
		      GTK_DEST_DEFAULT_ALL,
		      targ, 1,
		      GDK_ACTION_COPY);

    g_signal_connect(G_OBJECT(w), "drag-data-received",
		     G_CALLBACK(session_data_received),
		     NULL);
}

static void drag_graph (GtkWidget *w, GdkDragContext *context,
			GtkSelectionData *sel, guint info, guint t,
			SESSION_GRAPH *graph)
{
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_STRING, 8,
                           (const guchar *) graph->fname,
			   strlen(graph->fname));
}

static void graph_drag_connect (GtkWidget *w, SESSION_GRAPH *graph)
{
    gtk_drag_source_set(w, GDK_BUTTON1_MASK,
                        &gretl_drag_targets[GRETL_GRAPH_FILE],
                        1, GDK_ACTION_COPY);
    g_signal_connect(G_OBJECT(w), "drag-data-get",
                     G_CALLBACK(drag_graph), graph);
}

static void drag_model (GtkWidget *w, GdkDragContext *context,
			GtkSelectionData *sel, guint info, guint t,
			SESSION_MODEL *mod)
{
    drag_model_src = mod->ptr;
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_STRING, 8,
                           (const guchar *) &mod->ptr,
			   sizeof mod->ptr);
}

static void model_drag_connect (GtkWidget *w, SESSION_MODEL *mod)
{
    gtk_drag_source_set(w, GDK_BUTTON1_MASK,
                        &gretl_drag_targets[GRETL_MODEL_PTR],
                        1, GDK_ACTION_COPY);
    g_signal_connect(G_OBJECT(w), "drag-data-get",
                     G_CALLBACK(drag_model), mod);
}

#define WANT_TOOLTIP(t) (t == GRETL_OBJ_EQN || \
	                 t == GRETL_OBJ_GRAPH || \
                         t == GRETL_OBJ_PLOT)

static gui_obj *session_add_icon (gpointer data, int sort, int mode)
{
    gui_obj *obj;
    gchar *name = NULL;
    SESSION_GRAPH *graph = NULL;
    SESSION_MODEL *mod = NULL;
    SESSION_TEXT *text = NULL;
    int icon_named = 0;

    switch (sort) {
    case GRETL_OBJ_EQN:
    case GRETL_OBJ_VAR:
    case GRETL_OBJ_SYS:
	mod = (SESSION_MODEL *) data;
	name = g_strdup(mod->name);
	break;
    case GRETL_OBJ_PLOT:
    case GRETL_OBJ_GRAPH:
	graph = (SESSION_GRAPH *) data;
	name = g_strdup(graph->name);
	break;
    case GRETL_OBJ_TEXT:
	text = (SESSION_TEXT *) data;
	name = g_strdup(text->name);
	break;
    case GRETL_OBJ_DSET:
	name = g_strdup(_("Data set"));
	break;
    case GRETL_OBJ_SCALARS:
	name = g_strdup(_("Scalars"));
	break;
    case GRETL_OBJ_INFO:
	name = g_strdup(_("Data info"));
	break;
    case GRETL_OBJ_NOTES:
	name = g_strdup(_("Notes"));
	break;
    case GRETL_OBJ_CORR:
	name = g_strdup(_("Correlations"));
	break;
    case GRETL_OBJ_STATS:
	name = g_strdup(_("Summary"));
	break;
    case GRETL_OBJ_MODTAB:
	name = g_strdup(_("Model table"));
	break;
    case GRETL_OBJ_GPAGE:
	name = g_strdup(_("Graph page"));
	break;
    case GRETL_OBJ_MATRIX:
    case GRETL_OBJ_BUNDLE:
	name = g_strdup(user_var_get_name((user_var *) data));
	if (name == NULL || *name == '\0') {
	    fprintf(stderr, "session_add_icon: got no name for %s at %p\n",
		    sort == GRETL_OBJ_MATRIX ? "matrix" : "bundle",
		    (void *) data);
	}
	break;
    default:
	break;
    }

    if (name == NULL) {
	fprintf(stderr, "session_add_icon: NULL name for object of sort %d\n",
		sort);
	return NULL;
    }

    obj = gui_object_new(name, sort, data);

    /* full-length object name as tooltip */
    if (g_utf8_strlen(name, -1) > SHOWNAMELEN) {
	gretl_tooltips_add(GTK_WIDGET(obj->icon), name);
	icon_named = 1;
    }

    /* set up for drag and drop */
    if (sort == GRETL_OBJ_EQN) {
	model_drag_connect(obj->icon, obj->data);
    } else if (sort == GRETL_OBJ_GRAPH) {
	graph_drag_connect(obj->icon, obj->data);
    }

    /* second try at adding a tooltip */
    if (WANT_TOOLTIP(sort) && !icon_named) {
	char *str = NULL;

	if (sort == GRETL_OBJ_EQN) {
	    MODEL *pmod = (MODEL *) mod->ptr;

	    str = model_cmd_str(pmod);
	} else if (sort == GRETL_OBJ_GRAPH ||
		   sort == GRETL_OBJ_PLOT) {
	    str = graph_str(graph);
	}
	if (str != NULL) {
	    gretl_tooltips_add(GTK_WIDGET(obj->icon), str);
	    free(str);
	}
    }

    if (sort == GRETL_OBJ_DSET) {
	obj->data = datafile;
    }

    if (mode == ICON_ADD_SINGLE) {
	pack_single_icon(obj);
    } else if (mode == ICON_ADD_BATCH) {
	iconlist = g_list_append(iconlist, obj);
    }

    return obj;
}

static GtkWidget *create_pop_item (GtkWidget *popup, char *str,
				   GtkCallback callback)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(str);
    g_signal_connect(G_OBJECT(item), "activate",
		     G_CALLBACK(callback),
		     str);
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);

    if (!strcmp(str, _("Save session"))) {
	save_item = item;
	if (session.status == 0) {
	    gtk_widget_set_sensitive(item, FALSE);
	}
    }

    return item;
}

static void session_build_popups (void)
{
    GtkWidget *item;
    size_t i, n;

    mtab_add_item = NULL;
    gpage_add_item = NULL;
    graph_copy_item = NULL;

    if (global_popup == NULL) {
	global_popup = gtk_menu_new();
	n = G_N_ELEMENTS(global_items);
	for (i=0; i<n; i++) {
	    create_pop_item(global_popup, _(global_items[i]),
			    global_popup_callback);
	}
    }

    if (model_popup == NULL) {
	model_popup = gtk_menu_new();
	n = G_N_ELEMENTS(model_items);
	for (i=0; i<n; i++) {
	    item = create_pop_item(model_popup, _(model_items[i]),
				   object_popup_callback);
	    if (i == ADD_TO_MTAB_IDX) {
		mtab_add_item = item;
	    }
	}
    }

    if (model_table_popup == NULL) {
	model_table_popup = gtk_menu_new();
	n = G_N_ELEMENTS(model_table_items);
	for (i=0; i<n; i++) {
	    create_pop_item(model_table_popup, _(model_table_items[i]),
			    object_popup_callback);
	}
    }

    if (generic_popup == NULL) {
	generic_popup = gtk_menu_new();
	n = G_N_ELEMENTS(generic_items);
	for (i=0; i<n; i++) {
	    create_pop_item(generic_popup, _(generic_items[i]),
			    object_popup_callback);
	}
    }

    if (graph_popup == NULL) {
	graph_popup = gtk_menu_new();
	n = G_N_ELEMENTS(graph_items);
	for (i=0; i<n; i++) {
	    item = create_pop_item(graph_popup, _(graph_items[i]),
				   object_popup_callback);
	    if (i == ADD_TO_GPAGE_IDX) {
		gpage_add_item = item;
	    } else if (i == GRAPH_COPY_IDX) {
		graph_copy_item = item;
	    }
	}
    }

    if (graph_page_popup == NULL) {
	graph_page_popup = gtk_menu_new();
	n = G_N_ELEMENTS(graph_page_items);
	for (i=0; i<n; i++) {
	    create_pop_item(graph_page_popup, _(graph_page_items[i]),
			    object_popup_callback);
	}
    }

    if (data_popup == NULL) {
	data_popup = gtk_menu_new();
	n = G_N_ELEMENTS(dataset_items);
	for (i=0; i<n; i++) {
	    create_pop_item(data_popup, _(dataset_items[i]),
			    data_popup_callback);
	}
    }

    if (scalars_popup == NULL) {
	scalars_popup = gtk_menu_new();
	n = G_N_ELEMENTS(scalars_items);
	for (i=0; i<n; i++) {
	    create_pop_item(scalars_popup, _(scalars_items[i]),
			    scalars_popup_callback);
	}
    }

    if (info_popup == NULL) {
	info_popup = gtk_menu_new();
	n = G_N_ELEMENTS(info_items);
	for (i=0; i<n; i++) {
	    create_pop_item(info_popup, _(info_items[i]),
			    info_popup_callback);
	}
    }

    if (matrix_popup == NULL) {
	matrix_popup = gtk_menu_new();
	n = G_N_ELEMENTS(matrix_items);
	for (i=0; i<n; i++) {
	    create_pop_item(matrix_popup, _(matrix_items[i]),
			    matrix_popup_callback);
	}
    }

    if (bundle_popup == NULL) {
	bundle_popup = gtk_menu_new();
	n = G_N_ELEMENTS(bundle_items);
	for (i=0; i<n; i++) {
	    create_pop_item(bundle_popup, _(bundle_items[i]),
			    bundle_popup_callback);
	}
    }
}

#if 0

/* The reorganization of icons carried out by the following callback
   is not at all intuitive -- 2022-07-12
*/

static gboolean
iconview_resize_callback (GtkWidget *w, GdkEventConfigure *e, gpointer p)
{
    if (e->width != iconview_width) {
	if (iconview_width > 0) {
	    int cols = e->width / 100;

	    if (cols >= ICONVIEW_MIN_COLS && cols != iconview_cols) {
		iconview_cols = cols;
		rearrange_icons();
	    }
	}
	iconview_width = e->width;
    }

    return FALSE;
}

#endif

void view_session (void)
{
    GtkWidget *ebox, *scroller;
    gchar *title;
    int hmax = get_screen_height() / 2;
    int hmin = 360;
    int height;

    if (iconview != NULL) {
	gtk_window_present(GTK_WINDOW(iconview));
	return;
    }

    session_view_init();

    iconview = gretl_gtk_window();
    gtk_window_set_position(GTK_WINDOW(iconview), GTK_WIN_POS_MOUSE);
    title = g_strdup_printf("gretl: %s", _("icon view"));
    gtk_window_set_title(GTK_WINDOW(iconview), title);
    g_free(title);

    gtk_container_set_border_width(GTK_CONTAINER(iconview), 0);
    g_signal_connect(G_OBJECT(iconview), "destroy",
		     G_CALLBACK(session_view_free), NULL);

    session_build_popups();

    ebox = gtk_event_box_new();
    gtk_container_set_border_width(GTK_CONTAINER(ebox), 5);
    gtk_container_add(GTK_CONTAINER(iconview), ebox);
    g_signal_connect(G_OBJECT(ebox), "button-press-event",
		     G_CALLBACK(session_view_click), NULL);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_container_set_border_width(GTK_CONTAINER(scroller), 0);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(ebox), scroller);

    icon_table = gtk_table_new(2, ICONVIEW_MIN_COLS, FALSE);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller),
					  icon_table);

    height = 90 * add_all_icons();
    if (height < hmin) {
	height = hmin;
    } else if (height > hmax) {
	height = hmax;
    }
    gtk_window_set_default_size(GTK_WINDOW(iconview), 640, height);

    window_list_add(iconview, OPEN_SESSION);
    g_signal_connect(G_OBJECT(iconview), "key-press-event",
		     G_CALLBACK(catch_iconview_key), NULL);
#if 0 /* scrubbed 2022-07-12 */
    g_signal_connect(G_OBJECT(iconview), "configure-event",
		     G_CALLBACK(iconview_resize_callback), NULL);
#endif

    gtk_widget_show_all(iconview);

    gtk_container_foreach(GTK_CONTAINER(scroller),
			  (GtkCallback) white_bg_style,
			  NULL);

    gtk_widget_set_can_focus(icon_table, TRUE);
    gtk_widget_grab_focus(icon_table);
}

static int view_session_deferred;

static void auto_view_session (void)
{
    if (waiting_for_output()) {
	view_session_deferred = 1;
    } else {
	view_session();
    }
}

void maybe_view_session (void)
{
    if (view_session_deferred) {
	view_session_deferred = 0;
	view_session();
    }
}

static void make_short_label_string (char *targ, const char *src)
{
    if (g_utf8_strlen(src, -1) > SHOWNAMELEN) {
	g_utf8_strncpy(targ, src, SHOWNAMELEN - 3);
	strcat(targ, "...");
    } else {
	strcpy(targ, src);
    }
}

#ifdef REMEDY_LABELS

static gchar *rewrite_icon_label (const char *s)
{
    gchar *ret = NULL;
    int i, l1, l2, ld, ldmin = 100;
    int pos = 0;

    for (i=0; s[i] != '\0'; i++) {
	if (s[i] == ' ') {
	    l1 = g_utf8_strlen(s, i);
	    l2 = g_utf8_strlen(s + i + 1, -1);
	    ld = abs(l1 - l2);
	    if (ld < ldmin) {
		ldmin = ld;
		pos = i;
	    }
	}
    }

    if (pos > 0) {
	ret = g_strdup(s);
	ret[pos] = '\n';
    }

    return ret;
}

#endif

static void create_gobj_icon (gui_obj *obj, const char **xpm)
{
    GdkPixbuf *pbuf;
    GtkWidget *image;

    pbuf = gdk_pixbuf_new_from_xpm_data(xpm);

    obj->icon = gtk_event_box_new();
    gtk_widget_set_size_request(obj->icon, 44, 36);
    obj->label = NULL;

    image = gtk_image_new_from_pixbuf(pbuf);
    g_object_unref(G_OBJECT(pbuf));

    gtk_container_add(GTK_CONTAINER(obj->icon), image);
    gtk_widget_show(image);

    if (obj->sort == GRETL_OBJ_MODTAB || obj->sort == GRETL_OBJ_GPAGE) {
	session_drag_setup(obj);
    }

#ifdef REMEDY_LABELS
    if (g_utf8_strlen(obj->name, -1) > 12 && strchr(obj->name, ' ')) {
	gchar *tmp = rewrite_icon_label(obj->name);

	if (tmp != NULL) {
	    obj->label = gtk_label_new(tmp);
	    g_free(tmp);
	}
    }
#endif
    if (obj->label == NULL) {
	obj->label = gtk_label_new(obj->name);
    }
    gtk_label_set_width_chars(GTK_LABEL(obj->label), 12);
    gtk_label_set_max_width_chars(GTK_LABEL(obj->label), SHOWNAMELEN);
    gtk_label_set_line_wrap(GTK_LABEL(obj->label), TRUE);
    gtk_label_set_justify(GTK_LABEL(obj->label), GTK_JUSTIFY_CENTER);

    g_object_ref(obj->icon);
    g_object_ref(obj->label);

    g_signal_connect(G_OBJECT(obj->icon), "button-press-event",
		     G_CALLBACK(session_icon_click), obj);
    g_signal_connect(G_OBJECT(obj->icon), "enter-notify-event",
		     G_CALLBACK(icon_entered), obj);
    g_signal_connect(G_OBJECT(obj->icon), "leave-notify-event",
		     G_CALLBACK(icon_left), obj);
}

static gui_obj *gui_object_new (gchar *name, int sort, gpointer data)
{
    gui_obj *obj;
    char **xpm = NULL;

    obj = mymalloc(sizeof *obj);
    if (obj == NULL) {
	return NULL;
    }

    obj->name = name;
    obj->sort = sort;
    obj->data = data;

#if SESSION_DEBUG
    fprintf(stderr, "Allocated obj at %p (%s)\n", (void *) obj, obj->name);
#endif

    switch (sort) {
    case GRETL_OBJ_EQN:
    case GRETL_OBJ_VAR:
    case GRETL_OBJ_SYS:     xpm = model_xpm;       break;
    case GRETL_OBJ_PLOT:    xpm = boxplot_xpm;     break;
    case GRETL_OBJ_GRAPH:   xpm = gnuplot_xpm;     break;
    case GRETL_OBJ_DSET:    xpm = dot_sc_xpm;      break;
    case GRETL_OBJ_SCALARS: xpm = dot_sc_xpm;      break;
    case GRETL_OBJ_INFO:    xpm = xfm_info_xpm;    break;
    case GRETL_OBJ_TEXT:
    case GRETL_OBJ_NOTES:   xpm = text_xpm;        break;
    case GRETL_OBJ_CORR:    xpm = rhohat_xpm;      break;
    case GRETL_OBJ_STATS:   xpm = summary_xpm;     break;
    case GRETL_OBJ_MODTAB:  xpm = model_table_xpm; break;
    case GRETL_OBJ_GPAGE:   xpm = graph_page_xpm;  break;
    case GRETL_OBJ_MATRIX:  xpm = matrix_xpm;      break;
    case GRETL_OBJ_BUNDLE:  xpm = bundle_xpm;      break;
    default: break;
    }

    create_gobj_icon(obj, (const char **) xpm);

    return obj;
}

static void real_open_session_graph (SESSION_GRAPH *graph)
{
    GtkWidget *plotwin = get_window_for_plot(graph);

#if GRAPH_DEBUG
    fprintf(stderr, "real_open_session_graph: graph %p, plotwin %p\n",
	    (void *) graph, (void *) plotwin);
#endif

    if (plotwin != NULL) {
	gtk_window_present(GTK_WINDOW(plotwin));
    } else {
	char tmp[MAXLEN];

	session_file_make_path(tmp, graph->fname, NULL);
#if GRAPH_DEBUG
	fprintf(stderr, " graph->fname = '%s'\n", graph->fname);
	fprintf(stderr, " tmp = '%s'\n", tmp);
	fprintf(stderr, " call display_session_graph\n");
#endif
	display_session_graph(tmp, graph->name, graph);
    }
}

static void open_gui_graph (gui_obj *obj)
{
    real_open_session_graph((SESSION_GRAPH *) obj->data);
}

void display_session_graph_by_data (void *p)
{
    real_open_session_graph((SESSION_GRAPH *) p);
}

gchar *session_graph_get_filename (void *p)
{
    if (p != NULL) {
	SESSION_GRAPH *graph = p;
	char tmp[MAXLEN];

	session_file_make_path(tmp, graph->fname, NULL);
	return g_strdup(tmp);
    } else {
	return NULL;
    }
}

static int is_idempotent (const gretl_matrix *m,
			  const gretl_matrix *evals)
{
    double tol = 1.0e-12;

    if (evals != NULL) {
	double x;
	int i;

	if (gretl_is_complex(evals)) {
	    double complex z;

	    for (i=0; i<m->rows; i++) {
		z = gretl_cmatrix_get(evals, i, 0);
		if (fabs(cimag(z)) > tol) {
		    return 0;
		} else {
		    x = creal(z);
		    if (fabs(x * (1.0 - x)) > tol) {
			return 0;
		    }
		}
	    }
	} else {
	    for (i=0; i<m->rows; i++) {
		x = evals->val[i];
		if (fabs(x * (1.0 - x)) > tol) {
		    return 0;
		}
	    }
	}
    }

    return gretl_matrix_is_idempotent(m, tol);
}

static void print_int_formatted (char *s, int k, PRN *prn)
{
    int len = 12, n = strlen(s) - g_utf8_strlen(s, -1);
    char fmt[24];

    if (n > 0) {
	len += n;
    }

    sprintf(fmt, "%%-%ds %%3d\n", len);
    pprintf(prn, fmt, s, k);
}

static void print_double_formatted (char *s, double x, PRN *prn)
{
    int len = 16, n = strlen(s) - g_utf8_strlen(s, -1);
    char fmt[24];

    if (n > 0) {
	len += n;
    }

    sprintf(fmt, "%%-%ds %%.8g\n", len);
    pprintf(prn, fmt, s, x);
}

void
view_matrix_properties (const gretl_matrix *m, const char *name)
{
    gretl_matrix *A = NULL;
    gretl_matrix *evals = NULL;
    gchar *title;
    PRN *prn;
    int s, err = 0;

    if (m == NULL || bufopen(&prn)) {
	return;
    }

    pprintf(prn, _("Properties of matrix %s"), (name != NULL)? name : "");
    pputs(prn, "\n\n");

    if (m->rows == 0 || m->cols == 0) {
	pprintf(prn, _("Null matrix, %d x %d\n"), m->rows, m->cols);
	goto done;
    } else if (m->rows == 1 && m->cols == 1) {
	pprintf(prn, _("Scalar matrix, value %g\n"), m->val[0]);
	goto done;
    } else if (gretl_is_identity_matrix(m)) {
	pprintf(prn, _("Identity matrix, order %d\n"), m->rows);
	goto done;
    } else if (gretl_is_zero_matrix(m)) {
	pprintf(prn, _("Null matrix, %d x %d\n"), m->rows, m->cols);
	goto done;
    }

    print_int_formatted(_("Rows"), m->rows, prn);
    print_int_formatted(_("Columns"), m->cols, prn);
    print_int_formatted(_("Rank"), gretl_matrix_rank(m, NADBL, &err), prn);

    s = gretl_matrix_get_structure(m);

    if (s > 0) {
	pprintf(prn, "%s\n", _("Square"));
    }

    if (s == GRETL_MATRIX_DIAGONAL) {
	pprintf(prn, "%s\n", _("Diagonal"));
    } else if (s == GRETL_MATRIX_LOWER_TRIANGULAR) {
	pprintf(prn, "%s\n", _("Lower triangular"));
    } else if (s == GRETL_MATRIX_UPPER_TRIANGULAR) {
	pprintf(prn, "%s\n", _("Upper triangular"));
    } else if (s == GRETL_MATRIX_SYMMETRIC) {
	pprintf(prn, "%s\n", _("Symmetric"));
	A = gretl_matrix_copy(m);
	if (A != NULL) {
	    err = gretl_matrix_cholesky_decomp(A);
	    if (!err) {
		pprintf(prn, "%s\n", _("Positive definite"));
	    } else {
		pprintf(prn, "%s\n", _("Not positive definite"));
	    }
	    gretl_matrix_copy_values(A, m);
	    err = 0;
	    evals = gretl_symmetric_matrix_eigenvals(A, 0, &err);
	}
    }

    if (s > 0 && (s != GRETL_MATRIX_SYMMETRIC || evals == NULL)) {
	evals = gretl_dgeev(m, NULL, NULL, &err);
    }

    if (s > 0) {
	if (is_idempotent(m, evals)) {
	    pprintf(prn, "%s\n", _("Idempotent"));
	} else {
	    pprintf(prn, "%s\n", _("Not idempotent"));
	}
    }

    pputc(prn, '\n');

    print_double_formatted(_("1-norm"), gretl_matrix_one_norm(m), prn);
    print_double_formatted(_("Infinity-norm"), gretl_matrix_infinity_norm(m), prn);

    if (m->rows == m->cols) {
	double det;

	print_double_formatted(_("Trace"), gretl_matrix_trace(m), prn);
	if (A == NULL) {
	    A = gretl_matrix_copy(m);
	} else {
	    gretl_matrix_copy_values(A, m);
	}
	if (A != NULL) {
	    det = gretl_matrix_determinant(A, &err);
	    if (!err) {
		print_double_formatted(_("Determinant"), det, prn);
	    }
	}
    }

    if (evals != NULL) {
	double complex z;
	int i;

	pprintf(prn, "\n%s:\n", _("Eigenvalues"));

	for (i=0; i<m->rows; i++) {
	    if (gretl_is_complex(evals)) {
		z = gretl_cmatrix_get(evals, i, 0);
		pprintf(prn, "  (%.8g, %.8g)\n", creal(z), cimag(z));
	    } else {
		pprintf(prn, "  %.8g\n", evals->val[i]);
	    }
	}
	gretl_matrix_free(evals);
    }

    if (A != NULL) {
	gretl_matrix_free(A);
    }

 done:

    title = gretl_window_title(name);
    view_buffer(prn, 78, 400, title, PRINT, NULL);
    g_free(title);
}
