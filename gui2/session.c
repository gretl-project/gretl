/*
 *  Copyright (c) by Allin Cottrell
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

/* session.c for gretl */

#include "gretl.h"
#include "session.h"
#include "selector.h"
#include "boxplots.h"
#include "ssheet.h"
#include "gpt_control.h"
#include "guiprint.h"
#include "model_table.h"
#include "graph_page.h"
#include "textbuf.h"
#include "cmdstack.h"
#include "filelists.h"
#include "dlgutils.h"
#include "fileselect.h"

#include "var.h"
#include "varprint.h"
#include "objstack.h"
#include "system.h"

#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#ifdef _WIN32
# include <windows.h>
#endif

#ifdef OLD_GTK
# include <gdk-pixbuf/gdk-pixbuf.h>
#endif

static void auto_save_gp (windata_t *vwin);

#include "../pixmaps/model.xpm"
#include "../pixmaps/boxplot.xpm"
#include "../pixmaps/gnuplot.xpm"
#include "../pixmaps/xfm_sc.xpm"
#include "../pixmaps/xfm_info.xpm"
#include "../pixmaps/xfm_text.xpm"
#include "../pixmaps/xfm_make.xpm"
#include "../pixmaps/rhohat.xpm"
#include "../pixmaps/summary.xpm"
#include "../pixmaps/model_table.xpm"
#include "../pixmaps/graph_page.xpm"

#undef SESSION_DEBUG

#define OBJNAMLEN   32
#define SHOWNAMELEN 12

typedef struct SESSION_ SESSION;
typedef struct SESSIONBUILD_ SESSIONBUILD;
typedef struct SESSION_TEXT_ SESSION_TEXT;
typedef struct SESSION_MODEL_ SESSION_MODEL;
typedef struct gui_obj_ gui_obj;

struct SESSION_ {
    char name[MAXLEN];
    int nmodels;
    int ngraphs;
    int ntexts;
    GRAPHT **graphs;
    SESSION_MODEL **models;
    SESSION_TEXT **texts;
    char *notes;
};

struct SESSIONBUILD_ {
    int nmodels;
    int *model_ID;
    char **model_name;
};

struct SESSION_TEXT_ {
    char name[OBJNAMLEN];
    char *buf;
};

struct SESSION_MODEL_ {
    char *name;
    void *ptr;
    GretlObjType type;
};

struct gui_obj_ {
    gchar *name;
    gint sort;
    gpointer data;
    GtkWidget *icon;
    GtkWidget *label;
    gint row, col;
};

enum {
    ICON_ADD_BATCH,
    ICON_ADD_SINGLE
} icon_add_modes;

enum {
    SAVEFILE_SESSION,
    SAVEFILE_SCRIPT,
    SAVEFILE_ERROR
} savefile_returns;

static GtkTargetEntry session_drag_targets[] = {
    { "model_pointer", GTK_TARGET_SAME_APP, GRETL_MODEL_POINTER },
    { "text/uri-list", 0, GRETL_FILENAME }
};

static char *global_items[] = {
    N_("Arrange icons"),
    N_("Close window")
};

static char *model_items[] = {
    N_("Display"),
    N_("Add to model table"),
    N_("Delete")
};

static char *model_table_items[] = {
    N_("Display"),
    N_("Clear"),
    N_("Options"),
    N_("Help")
};

static char *graph_page_items[] = {
    N_("Display"),
    N_("Save as TeX..."),
    N_("Clear"),
    N_("Help")
};

static char *var_items[] = {
    N_("Display"),
    N_("Delete")
};

static char *graph_items[] = {
    N_("Display"),
    N_("Edit plot commands"),
    N_("Delete")
};

static char *dataset_items[] = {
    N_("Edit"),
    N_("Save..."),
    N_("Export as CSV..."),
    N_("Copy as CSV...")
};

static char *info_items[] = {
    N_("View"),
    N_("Edit"),
};

static char *session_items[] = {
    N_("Save"),
    N_("Save As...")
};

/* file-scope globals */

SESSION session;            /* hold models, graphs */
SESSIONBUILD rebuild;       /* rebuild session later */

static int session_file_open;

static GtkWidget *iconview;
static GtkWidget *icon_table;
static GtkWidget *global_popup;
static GtkWidget *session_popup;
static GtkWidget *model_popup;
static GtkWidget *model_table_popup;
static GtkWidget *var_popup;
static GtkWidget *graph_popup;
static GtkWidget *graph_page_popup;
static GtkWidget *boxplot_popup;
static GtkWidget *data_popup;
static GtkWidget *info_popup;

static GList *icon_list;
static gui_obj *active_object;

/* private functions */
static gui_obj *gui_object_new (gchar *name, int sort);
static gui_obj *session_add_icon (gpointer data, int sort, int mode);
static void session_build_popups (void);
static void global_popup_activated (GtkWidget *widget, gpointer data);
static void session_popup_activated (GtkWidget *widget, gpointer data);
static void object_popup_activated (GtkWidget *widget, gpointer data);
static void data_popup_activated (GtkWidget *widget, gpointer data);
static void info_popup_activated (GtkWidget *widget, gpointer data);
static void session_delete_icon (gui_obj *gobj);
static void open_gui_graph (gui_obj *gobj);
static void open_boxplot (gui_obj *gobj);
static gboolean session_icon_click (GtkWidget *widget, 
				    GdkEventButton *event,
				    gpointer data);
static void free_session_text (SESSION_TEXT *text);
static void free_session_model (SESSION_MODEL *mod);
static int real_delete_model_from_session (SESSION_MODEL *model);
#ifndef OLD_GTK
static void rename_session_object (gui_obj *obj, const char *newname);
#endif

#ifdef SESSION_DEBUG
static void print_session (const char *msg)
{
    int i;

    if (msg != NULL) {
	fprintf(stderr, "%s:\n", msg);
    }

    fprintf(stderr, "Session contains %d models\n", session.nmodels);
    for (i=0; i<session.nmodels; i++) {
	fprintf(stderr, " model '%s'\n", session.models[i]->name);
    }
    fprintf(stderr, "Session contains %d graphs\n", session.ngraphs);
    for (i=0; i<session.ngraphs; i++) {
	fprintf(stderr, " graph: %s (%s)\n", session.graphs[i]->name,
		session.graphs[i]->fname);
    }
    fprintf(stderr, "Session contains %d texts\n", session.ntexts);
    for (i=0; i<session.ntexts; i++) {
	fprintf(stderr, " text: '%s'\n", session.texts[i]->name);
    }
}
#endif

#if 0
static char *print_session_xml (void)
{
    PRN *prn;
    char *buf;
    int i;

    /* open the prn */

    pputs(prn, "<gui-session>\n");
    pprintf(prn, " <models number=\"%d\">\n", session.nmodels);
    for (i=0; i<session.nmodels; i++) {
	pprintf(prn, "  <session-model name=\"%s\" addr=\"%d\"/>\n", 
		session.models[i]->name, session.models[i]->ptr);
    }
    pputs(prn, " </models>\n");
    pprintf(prn, " <graphs number=\"%d\">\n", session.ngraphs);
    for (i=0; i<session.ngraphs; i++) {
	pprintf(prn, "  <session-graph name=\"%s\" fname=\"%s\"/>\n", 
		session.graphs[i]->name, session.graphs[i]->fname);
    } 
    pputs(prn, " </graphs>\n");
    pprintf(prn, " <texts number=\"%d\">\n", session.ntexts);
    for (i=0; i<session.ntexts; i++) {
	pprintf(prn, "  <session-text name=\"%s\">\n", session.texts[i]->name);
	/* XML encoding? */
	pputs(prn, session.texts[i]->buf);
	pputs(prn, "  </session-text>\n");
    }    
    pputs(prn, " </texts>\n");
    pputs(prn, "</gui-session>\n");

    /* close prn and return buf */
}
#endif

static int session_saved;

int session_is_saved (void)
{
    return session_saved;
}

void set_session_saved (int val)
{
    session_saved = val;
}

static void rebuild_init (void)
{
    rebuild.nmodels = 0;

    rebuild.model_ID = NULL;
    rebuild.model_name = NULL;
}

static void free_rebuild (void)
{
    int i;

    if (rebuild.model_ID != NULL) {
	free(rebuild.model_ID);
	rebuild.model_ID = NULL;
    }

    if (rebuild.model_name != NULL) {
	for (i=0; i<rebuild.nmodels; i++) {
	    free(rebuild.model_name[i]);
	}
	free(rebuild.model_name);
	rebuild.model_name = NULL;
    }
}

static void edit_session_notes (void)
{
    static GtkWidget *notes_window;

    if (notes_window == NULL) {
	windata_t *vwin;

	vwin = edit_buffer(&session.notes, 80, 400, 
			   _("gretl: session notes"),
			   EDIT_NOTES);
	notes_window = vwin->dialog;
	g_signal_connect(G_OBJECT(vwin->dialog), "destroy",
			 G_CALLBACK(gtk_widget_destroyed),
			 &notes_window);
    } else {
	gdk_window_show(notes_window->window);
	gdk_window_raise(notes_window->window);
    }
}

static int look_up_graph_by_name (const char *grname)
{
    int i;

    for (i=0; i<session.ngraphs; i++) {
	if (!strcmp(grname, (session.graphs[i])->name)) {
	    return i;
	}
    }
    return -1;
}

static int look_up_text_by_name (const char *tname)
{
    int i;

    for (i=0; i<session.ntexts; i++) {
	if (!strcmp(tname, session.texts[i]->name)) {
	    return i;
	}
    }
    return -1;
}

int real_add_graph_to_session (const char *fname, const char *grname,
			       int code)
{
    int ng = look_up_graph_by_name(grname);
    int replace = 0;

    if (ng >= 0) {
	replace = 1;
	(session.graphs[ng])->sort = code;
	strcpy((session.graphs[ng])->fname, fname);	
	(session.graphs[ng])->ID = plot_count++;
    } else {
	GRAPHT **graphs;

	ng = session.ngraphs;

	graphs = myrealloc(session.graphs, (ng + 1) * sizeof *graphs);
	if (graphs == NULL) {
	    return ADD_OBJECT_FAIL;
	}

	session.graphs = graphs;

	session.graphs[ng] = mymalloc(sizeof **session.graphs);
	if (session.graphs[ng] == NULL) {
	    return ADD_OBJECT_FAIL;
	}

	(session.graphs[ng])->sort = code;

	strcpy((session.graphs[ng])->fname, fname);
	strcpy((session.graphs[ng])->name, grname);
	(session.graphs[ng])->ID = plot_count++;
	session.ngraphs += 1;
    }
    
    session_changed(1);

    if (icon_list != NULL && !replace) {
	session_add_icon(session.graphs[ng], 
			 (code == GRETL_GNUPLOT_GRAPH)? GRETL_OBJ_GRAPH : GRETL_OBJ_PLOT,
			 ICON_ADD_SINGLE);
    }

    return (replace)? ADD_OBJECT_REPLACE : ADD_OBJECT_OK;
}

void add_graph_to_session (gpointer data, guint code, GtkWidget *w)
{
    char pltname[MAXLEN], savedir[MAXLEN];
    char grname[OBJNAMLEN];
    int boxplot_count;
    
    errno = 0;

    get_default_dir(savedir, SAVE_THIS_GRAPH);

    if (code == GRETL_GNUPLOT_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;

	sprintf(pltname, "%ssession.Graph_%d", savedir, plot_count + 1);
	sprintf(grname, "%s %d", _("Graph"), plot_count + 1);
	/* move temporary plot file to permanent */
	if (copyfile(plot->fname, pltname)) {
	    return;
	} 
	if (remove_png_term_from_plotfile(pltname, plot)) {
	    errbox(_("Failed to copy graph file"));
	    return;
	}
	remove(plot->fname);
	strcpy(plot->fname, pltname);
	mark_plot_as_saved(plot);	
    } else if (code == GRETL_BOXPLOT) {
	boxplot_count = augment_boxplot_count();
	sprintf(pltname, "%ssession.Plot_%d", savedir, boxplot_count);
	sprintf(grname, "%s %d", _("Boxplot"), boxplot_count);
	if (copyfile(boxplottmp, pltname)) {
	    return;
	} 
	remove(boxplottmp);
    } else {
	errbox("bad code in add_graph_to_session");
	return;
    }

    real_add_graph_to_session(pltname, grname, code);
}

static SESSION_MODEL *
session_model_new (void *ptr, const char *name, GretlObjType type)
{
    SESSION_MODEL *mod = malloc(sizeof *mod);

    if (mod != NULL) {
	mod->ptr = ptr;
	mod->type = type;
	if (name == NULL) {
	    gretl_object_compose_name(ptr, type);
	    name = gretl_object_get_name(ptr, type);
	} 
	mod->name = g_strdup(name);
    }

    return mod;
}

static int real_add_model_to_session (void *ptr, const char *name,
				      GretlObjType type)
{
    SESSION_MODEL **models;
    SESSION_MODEL *newmod;
    int nm = session.nmodels; 

    newmod = session_model_new(ptr, name, type);
    if (newmod == NULL) {
	return 1;
    }    

    models = myrealloc(session.models, (nm + 1) * sizeof *models);
    if (models == NULL) {
	free_session_model(newmod);
	return 1;
    }

    session.models = models;
    session.models[nm] = newmod;
    session.nmodels += 1;

    /* note: augment refcount for this model */
    gretl_object_ref(ptr, type);

    /* add model icon to session display */
    if (icon_list != NULL) {
	session_add_icon(session.models[nm], type, ICON_ADD_SINGLE);
    }

    return 0;
}

int real_add_text_to_session (PRN *prn, const char *tname)
{
    const char *pbuf;
    int nt = look_up_text_by_name(tname);
    int replace = 0;

    if (nt >= 0) {
	pbuf = gretl_print_get_buffer(prn);
	free(session.texts[nt]->buf);
	session.texts[nt]->buf = g_strdup(pbuf);
	replace = 1;
    } else {
	SESSION_TEXT **texts;

	nt = session.ntexts;

	texts = myrealloc(session.texts, (nt + 1) * sizeof *texts);
	if (texts == NULL) {
	    return ADD_OBJECT_FAIL;
	}

	session.texts = texts;

	session.texts[nt] = mymalloc(sizeof **session.texts);
	if (session.texts[nt] == NULL) {
	    return ADD_OBJECT_FAIL;
	}

	pbuf = gretl_print_get_buffer(prn);
	session.texts[nt]->buf = g_strdup(pbuf);

	strcpy(session.texts[nt]->name, tname);

	session.ntexts += 1;
    }
    
    session_changed(1);

    if (icon_list != NULL && !replace) {
	session_add_icon(session.texts[nt], GRETL_OBJ_TEXT, ICON_ADD_SINGLE);
    }

    return (replace)? ADD_OBJECT_REPLACE : ADD_OBJECT_OK;
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

static SESSION_MODEL *get_session_model_by_data (void *ptr)
{
    int i;

    for (i=0; i<session.nmodels; i++) {
	if (ptr == session.models[i]->ptr) {
	    return session.models[i];
	}
    }

    return NULL;
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

int maybe_add_model_to_session (void *ptr, GretlObjType type)
{
    SESSION_MODEL *oldmod;
    const char *name;

    if (get_session_model_by_data(ptr)) {
	return 1;
    }

    name = gretl_object_get_name(ptr, type);

    /* check to see if there's already a model with
       the same name: if so, delete it
    */
    oldmod = get_session_model_by_name(name);
    if (oldmod != NULL) {
	real_delete_model_from_session(oldmod);
    }

    return real_add_model_to_session(ptr, name, type);
}

void model_add_as_icon (gpointer p, guint type, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    void *ptr = vwin->data;

    if (ptr == NULL) {
	return;
    }

    if (get_session_model_by_data(ptr)) {
	infobox(_("Model is already saved"));
	return;
    }

    if (real_add_model_to_session(ptr, NULL, type)) {
	return;
    }

    session_changed(1);
}

void model_add_as_icon_and_close (gpointer p, guint type, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;

    model_add_as_icon(p, type, w);

    if (!window_is_busy(vwin)) {
	gtk_widget_destroy(gtk_widget_get_toplevel(GTK_WIDGET(vwin->w)));
    } 
}

int session_changed (int set)
{
    static int has_changed;
    int orig;

    orig = has_changed;
    if (set >= 0) {
	has_changed = set;
    }

    return orig;
}

void session_init (void)
{
    session.models = NULL;
    session.graphs = NULL;
    session.texts = NULL;
    session.notes = NULL;
    session.nmodels = 0;
    session.ngraphs = 0;
    session.ntexts = 0;
    *session.name = '\0';

    session_changed(0);
    winstack_init();
    session_file_manager(CLEAR_DELFILES, NULL);
}

static int pick_up_session_notes (const char *fname)
{
    char notesfile[MAXLEN];
    struct stat buf;

    switch_ext(notesfile, fname, "Notes");

    if (stat(notesfile, &buf) == 0) {
	char noteline[MAXLEN];
	FILE *fp = gretl_fopen(notesfile, "r");

	if (fp == NULL) {
	    return 1;
	}

	/* read into buffer */
	session.notes = mymalloc(buf.st_size + 1);
	if (session.notes == NULL) {
	    fclose(fp);
	    return 1;
	}

	session.notes[0] = '\0';
	while (fgets(noteline, sizeof noteline, fp)) {
	    strcat(session.notes, noteline);
	} 

	fclose(fp);
    }

    return 0;
}

static void
session_name_from_session_file (char *sname, const char *fname)
{
    const char *p = strrchr(fname, SLASH);
    char *q;

    if (p != NULL) {
	strcpy(sname, p + 1);
    } else {
	strcpy(sname, fname);
    }

    q = strstr(sname, ".gretl");
    if (q != NULL) {
	*q = '\0';
    }
}

void do_open_session (GtkWidget *w, gpointer data)
{
    dialog_t *d = NULL;
    windata_t *fwin = NULL;
    FILE *fp;
    int status;

    if (data != NULL) {    
	if (w == NULL) {
	    /* not coming from edit_dialog */
	    fwin = (windata_t *) data;
	} else {
	    d = (dialog_t *) data;
	    fwin = (windata_t *) edit_dialog_get_data(d);
	}
    }

    fp = gretl_fopen(tryscript, "r");
    if (fp != NULL) {
	fclose(fp);
	strcpy(scriptfile, tryscript);
    } else {
	gchar *errbuf;

	errbuf = g_strdup_printf(_("Couldn't open %s\n"), tryscript);
	errbox(errbuf);
	g_free(errbuf);
	delete_from_filelist(FILE_LIST_SESSION, tryscript);
	delete_from_filelist(FILE_LIST_SCRIPT, tryscript);
	return;
    }

    clear_data();
    free_session();
    session_init();

    fprintf(stderr, I_("\nReading session file %s\n"), scriptfile);

    status = parse_savefile(scriptfile);

    if (status == SAVEFILE_ERROR) {
	return;
    } else if (status == SAVEFILE_SCRIPT) {
	do_open_script();
	return;
    }

    if (recreate_session(scriptfile)) { 
	return;
    }

#ifdef SESSION_DEBUG
    print_session("after recreate_session()");
#endif

    mkfilelist(FILE_LIST_SESSION, scriptfile);

    session_name_from_session_file(session.name, scriptfile);

    if (status == SAVEFILE_SESSION) {
	pick_up_session_notes(scriptfile);
    }

    /* trash the practice files window that launched the query? */
    if (fwin != NULL) {
	gtk_widget_destroy(fwin->w);   
    } 

    /* sync gui with session */
    session_file_open = 1;
    session_menu_state(TRUE);

#ifdef SESSION_DEBUG
    print_session("about to view_session()");
#endif

    view_session();
}

void verify_clear_data (void)
{
    if (dataset_locked()) {
	return;
    }

    if (!expert) {
        if (yes_no_dialog ("gretl",                      
			   _("Clearing the data set will end\n"
			     "your current session.  Continue?"), 
			   0) != GRETL_YES) {
            return;
	}
    }

    close_session();
}

static void free_session_model (SESSION_MODEL *mod)
{
    /* remove a reference to this model */
    gretl_object_unref(mod->ptr, mod->type);
    free(mod->name);
    free(mod);
}

void free_session (void)
{
    int i;

    if (session.models) {
	for (i=0; i<session.nmodels; i++) {
	    free_session_model(session.models[i]);
	}
	free(session.models);
	session.models = NULL;
    }

    if (session.graphs) {
	free(session.graphs);
	session.graphs = NULL;
    }

    if (session.texts) {
	for (i=0; i<session.ntexts; i++) {
	    free_session_text(session.texts[i]);
	}	
	free(session.texts);
	session.texts = NULL;
    }

    if (session.notes) {
	free(session.notes);
	session.notes = NULL;
    }

    session.nmodels = 0;
    session.ngraphs = 0;
    session.ntexts = 0;

    *session.name = '\0';
}

int highest_numbered_variable_in_session (void)
{
    MODEL *pmod;
    GRETL_VAR *var;
    gretl_equation_system *sys;
    int i, mvm, vmax = 0;

    if (session.models) {
	for (i=0; i<session.nmodels; i++) {
	    if (session.models[i]->type == GRETL_OBJ_EQN) {
		pmod = session.models[i]->ptr;
		if (pmod != NULL) {
		    mvm = highest_numbered_var_in_model(pmod, datainfo);
		    if (mvm > vmax) {
			vmax = mvm;
		    }
		}
	    } else if (session.models[i]->type == GRETL_OBJ_VAR) {
		var = session.models[i]->ptr;
		if (var != NULL) {
		    mvm = gretl_VAR_get_highest_variable(var, datainfo);
		    if (mvm > vmax) {
			vmax = mvm;
		    }
		}		
	    } else if (session.models[i]->type == GRETL_OBJ_SYS) {
		sys = session.models[i]->ptr;
		if (sys != NULL) {
		    mvm = highest_numbered_var_in_system(sys, datainfo);
		    if (mvm > vmax) {
			vmax = mvm;
		    }
		}
	    }
	}
    }

    return vmax;
}

int session_file_is_open (void)
{
    return session_file_open;
}

void close_session (void)
{
    clear_data(); 
    free_session();
    clear_model_table(NULL);
    clear_graph_page();

    session_menu_state(FALSE);
    session_file_open = 0;
    *scriptfile = '\0';

    if (iconview != NULL) {
	gtk_widget_destroy(iconview);
    }

    session_changed(0);
    set_session_saved(0);

    winstack_destroy();
    clear_selector();

    plot_count = 0;
    zero_boxplot_count();
}

/* see if there are saved objects from a previous session */

int saved_objects (const char *fname)
{
    FILE *fp;
    char saveline[MAXLEN];
    int saves = 0;

#ifdef SESSION_DEBUG
    fprintf(stderr, "saved_objects: checking '%s'\n", fname);
#endif    

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return -1;
    }

    while (fgets(saveline, sizeof saveline, fp)) {
	if (strncmp(saveline, "(* saved objects:", 17) == 0) {
	    saves = 1;
	    break;
	}
    }

    fclose(fp);

    return saves;
}

/* check session graph for name and filename, and check that
   the specified graph file exists and is readable 
*/

static int check_session_graph (const char *line, 
				char **name, char **fname)
{
    size_t len, lenmin = 24;
    const char *p;
    FILE *fp;

    p = strchr(line, '"');
    if (p == NULL) {
	return 1;
    }
    
    len = strcspn(++p, "\"");
    if (len + 1 < lenmin) {
	lenmin = len + 1;
    }

    *name = malloc(lenmin);
    **name = '\0';
    strncat(*name, p, lenmin - 1);

    p = strchr(p, '"'); 
    if (p == NULL) {
	free(*name);
	return 1;
    }

    *fname = g_strdup(++p);

    top_n_tail(*fname);

    fp = gretl_fopen(*fname, "r");

    if (fp == NULL) {
	gchar *msg = 
	    g_strdup_printf(_("Warning: couldn't open graph file %s"), *fname);

	errbox(msg);
	g_free(msg);		      
	free(*name);
	free(*fname);
	return 1;
    }

    fclose(fp);

    return 0;
}

static int allocate_model_rebuilder (int nmodels)
{
    int *IDs;
    char **names;

    IDs = myrealloc(rebuild.model_ID, nmodels * sizeof *IDs);
    if (IDs == NULL) {
	return E_ALLOC;
    }

    rebuild.model_ID = IDs;

    names = myrealloc(rebuild.model_name, nmodels * sizeof *names);
    if (names == NULL) {
	return E_ALLOC;
    }

    rebuild.model_name = names;

    names[nmodels - 1] = mymalloc(OBJNAMLEN);
    if (names[nmodels - 1] == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static int allocate_session_graph (int ngraphs)
{
    GRAPHT **graphs;

    graphs = myrealloc(session.graphs, ngraphs * sizeof *graphs);
    if (graphs == NULL) {
	return E_ALLOC;
    }

    session.graphs = graphs;

    graphs[ngraphs - 1] = mymalloc(sizeof **graphs);
    if (graphs[ngraphs - 1] == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static int saved_object_type (const char *str)
{
    if (!strcmp(str, "model")) {
	return GRETL_OBJ_EQN;
    } else if (!strcmp(str, "graph")) {
	return GRETL_OBJ_GRAPH;
    } else if (!strcmp(str, "plot")) {
	return GRETL_OBJ_PLOT;
    } else {
	return GRETL_OBJ_UNKNOWN;
    }
}

static int 
grab_model_from_session_file (const char *modline, int id)
{
    const char *p;
    int nm = rebuild.nmodels;
    int j, len;

#ifdef SESSION_DEBUG
    fprintf(stderr, "got a model to rebuild (%d)\n"
	    "rebuild.nmodels now -> %d\n", id, nm + 1);
#endif

    if (allocate_model_rebuilder(nm + 1)) {
	return 1;
    }

    rebuild.model_ID[nm] = id;

    p = strchr(modline, '"');
    if (p == NULL) {
	return 1;
    }

    p++;
    *rebuild.model_name[nm] = '\0';
    strncat(rebuild.model_name[nm], p, OBJNAMLEN - 1);
    len = strlen(rebuild.model_name[nm]);

    for (j=len; j>0; j--) {
	if (rebuild.model_name[nm][j] == '"') {
	    rebuild.model_name[nm][j] = '\0';
	    break;
	}
    }

#ifdef SESSION_DEBUG
    fprintf(stderr, "rebuild.model_name[%d] = '%s'\n", nm, 
	    rebuild.model_name[nm]);
#endif

    return 0;
}

static int 
grab_graph_from_session_file (const char *grline, int objtype)
{
    char *grname, *grfilename;
    int ng = session.ngraphs;

    if (check_session_graph(grline, &grname, &grfilename)) {
	return 1;
    }

    if (allocate_session_graph(ng + 1)) {
	return 1;
    }

    strcpy((session.graphs[ng])->name, grname);
    strcpy((session.graphs[ng])->fname, grfilename);

    free(grname);
    free(grfilename);

#ifdef SESSION_DEBUG
    fprintf(stderr, "got graph: '%s'\n", (session.graphs[ng])->fname);
#endif

    (session.graphs[ng])->ID = plot_count++;

    if (objtype == GRETL_OBJ_PLOT) {
	(session.graphs[ng])->sort = GRETL_BOXPLOT;
	augment_boxplot_count();
    } else {
	(session.graphs[ng])->sort = GRETL_GNUPLOT_GRAPH;
    }

    return 0;
}

int parse_savefile (const char *fname)
{
    FILE *fp;
    char fline[MAXLEN];
    int err = 0, objs = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return SAVEFILE_ERROR;
    }

    /* look for saved objects */
    while (fgets(fline, sizeof fline, fp)) {
	if (strncmp(fline, "(* saved objects:", 17) == 0) {
	    objs = 1;
	    break;
	}
    }

    if (objs == 0) { 
	/* no saved objects: just a regular script */
	fclose(fp);
	return SAVEFILE_SCRIPT;
    }

    rebuild_init();

#ifdef SESSION_DEBUG
    fprintf(stderr, "parse_savefile (%s): got saved objects\n", fname);
#endif

    while (fgets(fline, sizeof fline, fp) && !err) {
	int id, objtype = GRETL_OBJ_UNKNOWN;
	char objstr[8];

	if (strncmp(fline, "*)", 2) == 0) {
	    /* done reading saved objects */
	    break;
	}

	chopstr(fline);

	if (sscanf(fline, "%6s %d", objstr, &id) != 2) {
	    err = 1;
	} else {
	    objtype = saved_object_type(objstr);
	}

	/* FIXME types other than models, graphs? */

	if (objtype == GRETL_OBJ_EQN) {
	    if (grab_model_from_session_file(fline, id)) {
		err = 1;
	    } else {
		rebuild.nmodels += 1;
	    }
	} else if (objtype == GRETL_OBJ_GRAPH || objtype == GRETL_OBJ_PLOT) {
	    if (grab_graph_from_session_file(fline, objtype)) {
		err = 1;
	    } else {
		session.ngraphs += 1;
	    }
	} else {
	    err = 1;
	}
    }

    fclose(fp);

    if (err) {
	errbox(_("Session file is corrupted, ignoring"));
	return SAVEFILE_ERROR;
    }
    
    return SAVEFILE_SESSION;
}

/* called on start-up when a "session" file is loaded */

int recreate_session (const char *fname)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_NULL);

#ifdef SESSION_DEBUG
    fprintf(stderr, "recreate_session: fname = %s\n", fname);
    print_session(NULL);
#endif

    if (execute_script(fname, NULL, prn, REBUILD_EXEC)) {
	errbox(_("Error recreating session"));
    }

#ifdef SESSION_DEBUG
    print_session("recreate_session: after execute_script()");
#endif

    free_rebuild();
    gretl_print_destroy(prn);

    /* no fresh commands have been entered yet */
    set_replay_on();

    return 0;
}

void save_session_callback (GtkWidget *w, guint code, gpointer data)
{
    if (code == SAVE_AS_IS && session_file_open && scriptfile[0]) {
	save_session(scriptfile);
	session_changed(0);
    } else {
	file_selector(_("Save session"), SAVE_SESSION, FSEL_DATA_NONE, NULL);
    }
}

static char *model_cmd_str (MODEL *pmod)
{
    char *str = NULL;

    if (pmod->ci == MLE || pmod->ncoeff > 10 ||
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

static gchar *graph_str (GRAPHT *graph)
{
    FILE *fp;
    gchar *buf = NULL;

    fp = gretl_fopen(graph->fname, "r");

    if (fp != NULL) {
	char xlabel[24], ylabel[24], line[48];
	int gotxy = 0;

	while (fgets(line, 47, fp) && gotxy < 2) {
	    if (strstr(line, "# timeseries")) {
		break;
	    } else if (sscanf(line, "set xlabel %23s", xlabel) == 1) {
		gotxy++;
	    } else if (sscanf(line, "set ylabel %23s", ylabel) == 1) {
		gotxy++;
	    }
	}
	if (gotxy == 2) {
#ifdef OLD_GTK
	    buf = g_strdup_printf("%s %s %s", ylabel, _("versus"), xlabel);
#else
	    char *tmp = 
		g_strdup_printf("%s %s %s", ylabel, _("versus"), xlabel);

	    if (tmp != NULL) {
		buf = my_locale_to_utf8(tmp);
		free(tmp);
	    }
#endif /* OLD_GTK */
	}

	fclose(fp);
    }

    return buf;
}

static char *boxplot_str (GRAPHT *graph)
{
    FILE *fp;
    char *str = NULL;

    fp = gretl_fopen(graph->fname, "r");

    if (fp != NULL) {
	char vname[VNAMELEN], line[48];

	str = malloc(MAXLEN);
	if (str == NULL) {
	    return NULL;
	}
	*str = '\0';

	while (fgets(line, 47, fp) && strlen(str) < MAXLEN-48) {
	    chopstr(line);
	    if (sscanf(line, "%*d varname = %15s", vname) == 1) { 
		strcat(str, strchr(line, '=') + 2);
		strcat(str, " ");
	    }
	}
	fclose(fp);
    }

    return str;
}

static int maybe_raise_object_window (gpointer p)
{
    GtkWidget *w = match_window_by_data(p);
    int ret = 0;

    if (w != NULL) {
#ifdef OLD_GTK 
	gdk_window_show(w->window);
	gdk_window_raise(w->window);
#else
	gtk_window_present(GTK_WINDOW(w));
#endif
	ret = 1;
    }

    return ret;
}

static int display_session_model (SESSION_MODEL *sm)
{ 
    PRN *prn;
    int err = 0;

    if (maybe_raise_object_window(sm->ptr)) {
	return 0;
    }  

    if (bufopen(&prn)) {
	return 1;
    }    

    if (sm->type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) sm->ptr;

	printmodel(pmod, datainfo, OPT_NONE, prn);
	view_model(prn, pmod, 78, 400, sm->name);
    } else if (sm->type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) sm->ptr;

	gretl_VAR_print(var, datainfo, OPT_NONE, prn);
	view_buffer(prn, 78, 450, sm->name, var->ci, var);
    } else if (sm->type == GRETL_OBJ_SYS) {
	gretl_equation_system *sys = (gretl_equation_system *) sm->ptr;

	err = estimate_saved_equation_system(sys, &Z, datainfo, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    view_buffer(prn, 78, 450, sm->name, SYSTEM, sys);
	}
    }

    return err;
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
	/* FIXME not sure here: should we trash the icon too */
	gretl_object_unref(ptr, mod->type);
    } 
}	

static void open_boxplot (gui_obj *gobj)
{
    GRAPHT *graph = (GRAPHT *) gobj->data;

    retrieve_boxplot(graph->fname);
}

static void open_gui_text (gui_obj *gobj)
{ 
    SESSION_TEXT *text = (SESSION_TEXT *) gobj->data;
    PRN *prn;

    prn = gretl_print_new_with_buffer(g_strdup(text->buf));

    if (prn != NULL) { 
	view_buffer(prn, 80, 400, gobj->name, INFO, NULL);
    }
}

static void delete_delfiles (const gchar *fname, gpointer p)
{
    remove(fname);
}

void session_file_manager (int action, const char *fname)
{
    static GList *delfiles;

    if (action == SCHEDULE_FOR_DELETION) {
	gchar *delname;

	delname = g_strdup(fname);
	delfiles = g_list_append(delfiles, delname);
    } else if (action == REALLY_DELETE_ALL) {
	if (delfiles != NULL) {
	    g_list_foreach(delfiles, (GFunc) delete_delfiles, NULL);
	    g_list_free(delfiles);
	    delfiles = NULL;
	}
    } else if (action == CLEAR_DELFILES) {
	if (delfiles != NULL) {
	    g_list_free(delfiles);
	    delfiles = NULL;
	}
    }	     
}

static int real_delete_model_from_session (SESSION_MODEL *model)
{
    if (session.nmodels == 1) {
	free_session_model(session.models[0]);
    } else {
	SESSION_MODEL **mods;
	int i, j;

	mods = mymalloc((session.nmodels - 1) * sizeof *mods);
	if (session.nmodels > 1 && mods == NULL) {
	    return 1;
	}
	j = 0;
	for (i=0; i<session.nmodels; i++) {
	    if (session.models[i]->ptr != model->ptr) {
		mods[j++] = session.models[i];
	    } else {
		free_session_model(session.models[i]);
	    }
	}
	free(session.models);
	session.models = mods;
    }

    session.nmodels -= 1;
    session_changed(1);

    return 0;
}

static void free_session_text (SESSION_TEXT *text)
{
    free(text->buf);
    free(text);
}

static int real_delete_text_from_session (SESSION_TEXT *junk)
{
    if (session.ntexts == 1) {
	free_session_text(session.texts[0]);
    } else {
	SESSION_TEXT **pptext;
	int i, j;

	pptext = mymalloc((session.ntexts - 1) * sizeof *pptext);
	if (session.ntexts > 1 && pptext == NULL) {
	    return 1;
	}
	j = 0;
	for (i=0; i<session.ntexts; i++) {
	    if (session.texts[i] != junk) {
		pptext[j++] = session.texts[i];
	    } else {
		free_session_text(session.texts[i]);
	    }
	}
	free(session.texts);
	session.texts = pptext;
    }

    session.ntexts -= 1;
    session_changed(1);

    return 0;
}

static int real_delete_graph_from_session (GRAPHT *junk)
{
    int i, j;

    if (session.ngraphs == 1) {
	session_file_manager(SCHEDULE_FOR_DELETION,
			     (session.graphs[0])->fname);
	free(session.graphs[0]);
    } else {
	GRAPHT **ppgr;

	ppgr = mymalloc((session.ngraphs - 1) * sizeof *ppgr);
	if (session.ngraphs > 1 && ppgr == NULL) {
	    return 1;
	}
	j = 0;
	for (i=0; i<session.ngraphs; i++) {
	    if ((session.graphs[i])->ID != junk->ID) { 
		ppgr[j++] = session.graphs[i];
	    } else {
		session_file_manager(SCHEDULE_FOR_DELETION,
				     (session.graphs[i])->fname);
		free(session.graphs[i]);
	    }
	}
	free(session.graphs);
	session.graphs = ppgr;
    }

    session.ngraphs -= 1;
    session_changed(1);

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
    }

    set_replay_off();
    session_delete_icon(obj);

    return 0;
}

static void maybe_delete_session_object (gui_obj *obj)
{
    gpointer p;
    gchar *msg;

    if (obj->sort == GRETL_OBJ_EQN ||
	obj->sort == GRETL_OBJ_SYS || obj->sort == GRETL_OBJ_VAR) {
	SESSION_MODEL *mod = obj->data;

	p = mod->ptr;
    } else {
	p = obj->data;
    }

    if (winstack_match_data(p)) {
	errbox(_("Please close this object's window first"));
	return;
    }    

    msg = g_strdup_printf(_("Really delete %s?"), obj->name);
    if (yes_no_dialog(_("gretl: delete"), msg, 0) == GRETL_YES) {
	delete_session_object(obj);
    }

    g_free(msg);
}

#ifndef OLD_GTK

static void rename_session_graph (GRAPHT *graph, const char *newname)
{
    int i;

    for (i=0; i<session.ngraphs; i++) {
	if ((session.graphs[i])->ID == graph->ID) { 
	    (session.graphs[i])->name[0] = '\0';
	    strncat((session.graphs[i])->name, newname, 
		    sizeof (session.graphs[i])->name - 1);
	    break;
	}
    }
}

static void rename_session_object (gui_obj *obj, const char *newname)
{
    if (obj->sort == GRETL_OBJ_EQN || obj->sort == GRETL_OBJ_SYS || 
	obj->sort == GRETL_OBJ_VAR) { 
	SESSION_MODEL *sm = obj->data;

	gretl_object_rename(sm->ptr, sm->type, newname);
    } else if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) { 
	rename_session_graph(obj->data, newname);
    }

    free(obj->name);
    obj->name = g_strdup(newname);

    set_replay_off();
}

#endif /* !OLD_GTK */

static gui_obj *get_gui_obj_from_data (void *finddata)
{
    GList *mylist = icon_list;
    gui_obj *gobj = NULL;

    while (mylist != NULL) {
	gobj = (gui_obj *) mylist->data;
	if (gobj->data == finddata) {
	    return gobj;
	}
	mylist = mylist->next;
    }

    return NULL;
}

void delete_text_from_session (void *p)
{
    SESSION_TEXT *text = (SESSION_TEXT *) p;
    gui_obj *obj;

    if (text == NULL) return;

    real_delete_text_from_session(text);

    obj = get_gui_obj_from_data((void *) text);
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
    icon_list = NULL;
    icon_table = NULL;
}

static void delete_session_icon (gui_obj *gobj, gpointer p)
{
    if (gobj != NULL) {
	if (gobj->name != NULL) {
	    g_free(gobj->name); 
	}
	free(gobj);
    }
}

static void session_view_free (GtkWidget *w, gpointer data)
{
    iconview = NULL;

    g_list_foreach(icon_list, (GFunc) delete_session_icon, NULL);

    g_list_free(icon_list);
    icon_list = NULL;
}

#define SESSION_VIEW_COLS 4

static void session_delete_icon (gui_obj *gobj)
{
    if (gobj == NULL) return;

    if (gobj->icon && GTK_IS_WIDGET(gobj->icon))
	gtk_container_remove(GTK_CONTAINER(icon_table), gobj->icon);
    if (gobj->label && GTK_IS_WIDGET(gobj->label))
	gtk_container_remove(GTK_CONTAINER(icon_table), gobj->label);

    free(gobj->name);

    icon_list = g_list_first(icon_list);
    icon_list = g_list_remove(icon_list, gobj);
    if (icon_list == NULL) {
	fprintf(stderr, "Bad: icon_list has gone NULL\n");
    }
}

static void foreach_delete_icons (gui_obj *gobj, gpointer p)
{
    if (gobj->icon && GTK_IS_WIDGET(gobj->icon))
	gtk_container_remove(GTK_CONTAINER(icon_table), gobj->icon);
    if (gobj->label && GTK_IS_WIDGET(gobj->label))
	gtk_container_remove(GTK_CONTAINER(icon_table), gobj->label);

    free(gobj->name);
}

/* apparatus for getting a white background */

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

#ifdef OLD_GTK

static void white_bg_style (GtkWidget *widget, gpointer data)
{
     GtkRcStyle *rc_style;
     static GdkColor *white;

     if (white == NULL) {
	 white = get_white();
     }

     rc_style = gtk_rc_style_new();
     rc_style->bg[GTK_STATE_NORMAL] = *white;
     rc_style->color_flags[GTK_STATE_NORMAL] |= GTK_RC_BG;

     gtk_widget_modify_style(widget, rc_style);
     gtk_rc_style_unref(rc_style);
}

#else

static void white_bg_style (GtkWidget *widget, gpointer data)
{
    static GdkColor *white;

    if (white == NULL) {
	white = get_white();
    }

    gtk_widget_modify_bg(widget, GTK_STATE_NORMAL, white);
}

#endif /* OLD_GTK */


static void real_pack_icon (gui_obj *gobj, int row, int col)
{
    gobj->row = row;
    gobj->col = col;

    gtk_table_attach (GTK_TABLE(icon_table), gobj->icon,
		      col, col + 1, row, row + 1,
		      GTK_EXPAND, GTK_FILL, 5, 5);
    gtk_widget_show(gobj->icon);

    white_bg_style(gobj->icon, NULL);

    gtk_table_attach (GTK_TABLE(icon_table), gobj->label,
		      col, col + 1, row + 1, row + 2,
		      GTK_EXPAND, GTK_FILL, 5, 5);

    gtk_widget_show(gobj->label);
}

static void pack_single_icon (gui_obj *gobj)
{
    int row, col;
    gui_obj *last_obj;

    icon_list = g_list_last(icon_list);
    last_obj = icon_list->data;
    row = last_obj->row;
    col = last_obj->col;  

    icon_list = g_list_append(icon_list, gobj);

    col++;
    if (col > 0 && (col % SESSION_VIEW_COLS == 0)) {
	col = 0;
	row += 2;
	gtk_table_resize(GTK_TABLE(icon_table), 2 * row, SESSION_VIEW_COLS);
    } 

    real_pack_icon(gobj, row, col);
}

static void batch_pack_icons (void)
{
    gui_obj *gobj;
    int row = 0, col = 0;

    icon_list = g_list_first(icon_list);

    while (icon_list != NULL) {
	gobj = (gui_obj *) icon_list->data;
	real_pack_icon(gobj, row, col);
	col++;
	if (col > 0 && (col % SESSION_VIEW_COLS == 0)) {
	    col = 0;
	    row += 2;
	    gtk_table_resize(GTK_TABLE(icon_table), 2 * row, SESSION_VIEW_COLS);
	}
	if (icon_list->next == NULL) {
	    break;
	} else {
	    icon_list = icon_list->next;
	}
    }
}

static void add_all_icons (void) 
{
    int show_graph_page = check_for_prog(latex);
    int i;

    active_object = NULL;

    if (data_status) {
	session_add_icon(NULL, GRETL_OBJ_INFO,   ICON_ADD_BATCH);  /* data info */
	session_add_icon(NULL, GRETL_OBJ_DSET,   ICON_ADD_BATCH);  /* data file */
	session_add_icon(NULL, GRETL_OBJ_NOTES,  ICON_ADD_BATCH);  /* session notes */
	session_add_icon(NULL, GRETL_OBJ_STATS,  ICON_ADD_BATCH);  /* summary stats */
	session_add_icon(NULL, GRETL_OBJ_CORR,   ICON_ADD_BATCH);  /* correlation matrix */
	session_add_icon(NULL, GRETL_OBJ_MODTAB, ICON_ADD_BATCH);  /* model table */
	if (show_graph_page) {
	    session_add_icon(NULL, GRETL_OBJ_GPAGE, ICON_ADD_BATCH); /* graph page */
	}
    }

    session_add_icon(NULL, GRETL_OBJ_SCRIPT, ICON_ADD_BATCH);        /* script file */

#ifdef SESSION_DEBUG
    print_session("view_session");
#endif

    for (i=0; i<session.nmodels; i++) {
#ifdef SESSION_DEBUG
	fprintf(stderr, "adding session.models[%d] to view\n", i);
#endif
	session_add_icon(session.models[i], session.models[i]->type, 
			 ICON_ADD_BATCH);
    }

    for (i=0; i<session.ngraphs; i++) {
#ifdef SESSION_DEBUG
	fprintf(stderr, "adding session.graphs[%d] to view\n", i);
	fprintf(stderr, "this graph is of sort %d\n", (session.graphs[i])->sort);
#endif
	/* distinguish gnuplot graphs from gretl boxplots */
	session_add_icon(session.graphs[i], 
			 ((session.graphs[i])->sort == GRETL_BOXPLOT)? 
			 GRETL_OBJ_PLOT : GRETL_OBJ_GRAPH,
			 ICON_ADD_BATCH);
    }

    for (i=0; i<session.ntexts; i++) {
#ifdef SESSION_DEBUG
	fprintf(stderr, "adding session.texts[%d] to view\n", i);
#endif
	session_add_icon(session.texts[i], GRETL_OBJ_TEXT, ICON_ADD_BATCH);
    }

    batch_pack_icons();
}

static void rearrange_icons (void) 
{
    g_list_foreach(icon_list, (GFunc) foreach_delete_icons, NULL);

    g_list_free(icon_list);
    icon_list = NULL;

    add_all_icons();
}

static gint catch_iconview_key (GtkWidget *w, GdkEventKey *key, 
				gpointer p)
{
    if (key->keyval == GDK_q) { 
        gtk_widget_destroy(w);
    }

    return FALSE;
}

static void object_popup_show (gui_obj *gobj, GdkEventButton *event)
{
    GtkWidget *w = NULL;

    active_object = gobj;

    switch (gobj->sort) {
    case GRETL_OBJ_EQN: 
	w = model_popup; 
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
	w = var_popup; 
	break;
    case GRETL_OBJ_GRAPH: 
	w = graph_popup; 
	break;
    case GRETL_OBJ_PLOT: 
	w = boxplot_popup; 
	break;
    case GRETL_OBJ_DSET: 
	w = data_popup; 
	break;
    case GRETL_OBJ_INFO: 
	w = info_popup; 
	break;
    case GRETL_OBJ_SCRIPT: 
	w = session_popup; 
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
	errbox(_("The graph page is empty"));
    } else {
	file_selector(_("Save LaTeX file"), SAVE_TEX, FSEL_DATA_NONE, NULL);
    }
}

static void view_script_default (void)
{
    if (dump_command_stack(cmdfile, 0)) return;

    view_file(cmdfile, 0, 0, 78, 350, EDIT_SCRIPT);
}

static gboolean session_icon_click (GtkWidget *widget, 
				    GdkEventButton *event,
				    gpointer data)
{
    gui_obj *gobj;
    GdkModifierType mods;

    if (event == NULL) return FALSE;

    gdk_window_get_pointer(widget->window, NULL, NULL, &mods);

    if (data == NULL) { 
	/* click on window background */
	if (mods & GDK_BUTTON3_MASK) {
	    gtk_menu_popup(GTK_MENU(global_popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	}
	return TRUE;
    }

    gobj = (gui_obj *) data;

    if (event->type == GDK_2BUTTON_PRESS) {
	switch (gobj->sort) {
	case GRETL_OBJ_EQN:
	case GRETL_OBJ_VAR:
	case GRETL_OBJ_SYS:
	    display_session_model(gobj->data); break;
	case GRETL_OBJ_PLOT:
	    open_boxplot(gobj); break;
	case GRETL_OBJ_GRAPH:
	    open_gui_graph(gobj); break;
	case GRETL_OBJ_TEXT:
	    open_gui_text(gobj); break;
	case GRETL_OBJ_DSET:
	    show_spreadsheet(SHEET_EDIT_DATASET); break;
	case GRETL_OBJ_INFO:
	    open_info(NULL, 0, NULL); break;
	case GRETL_OBJ_SCRIPT:
	    view_script_default(); break;
	case GRETL_OBJ_NOTES:
	    edit_session_notes(); break;
	case GRETL_OBJ_MODTAB:
	    display_model_table_wrapper(); break;
	case GRETL_OBJ_GPAGE:
	    display_graph_page(); break;
	case GRETL_OBJ_CORR:
	    do_menu_op(NULL, CORR, NULL); break;
	case GRETL_OBJ_STATS:
	    do_menu_op(NULL, SUMMARY, NULL); break;
	}
	return TRUE;
    }

    if (mods & GDK_BUTTON3_MASK) {
	if (gobj->sort == GRETL_OBJ_EQN    || gobj->sort == GRETL_OBJ_GRAPH || 
	    gobj->sort == GRETL_OBJ_TEXT   || gobj->sort == GRETL_OBJ_DSET || 
	    gobj->sort == GRETL_OBJ_INFO   || gobj->sort == GRETL_OBJ_GPAGE ||
	    gobj->sort == GRETL_OBJ_SCRIPT || gobj->sort == GRETL_OBJ_PLOT || 
	    gobj->sort == GRETL_OBJ_MODTAB || gobj->sort == GRETL_OBJ_VAR ||
	    gobj->sort == GRETL_OBJ_SYS) {
	    object_popup_show(gobj, (GdkEventButton *) event);
	}
	return TRUE;
    }

    return FALSE;
}

static void global_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("Arrange icons")) == 0) 
	rearrange_icons();
    else if (strcmp(item, _("Close window")) == 0) 
	gtk_widget_destroy(iconview);
}

static void session_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("Save")) == 0) 
	save_session_callback(NULL, SAVE_AS_IS, NULL);
    else if (strcmp(item, _("Save As...")) == 0) 
	save_session_callback(NULL, SAVE_RENAME, NULL);
}

static void info_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("View")) == 0) 
	open_info(NULL, 0, NULL);
    else if (strcmp(item, _("Edit")) == 0) 
	edit_header(NULL, 0, NULL);
}

static void data_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("Edit")) == 0) 
	show_spreadsheet(SHEET_EDIT_DATASET);
    else if (strcmp(item, _("Save...")) == 0) 
	file_save(mdata, SAVE_DATA, NULL);
    else if (strcmp(item, _("Export as CSV...")) == 0) 
	file_save(mdata, EXPORT_CSV, NULL);
    else if (strcmp(item, _("Copy as CSV...")) == 0) 
	csv_to_clipboard();
}

static void object_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;
    gui_obj *obj;

    obj = active_object;

    if (strcmp(item, _("Display")) == 0) {
	if (obj->sort == GRETL_OBJ_EQN ||
	    obj->sort == GRETL_OBJ_VAR ||
	    obj->sort == GRETL_OBJ_SYS) {
	    display_session_model(obj->data);
	} else if (obj->sort == GRETL_OBJ_TEXT) {
	    open_gui_text(obj);
	} else if (obj->sort == GRETL_OBJ_MODTAB) {
	    display_model_table_wrapper();
	} else if (obj->sort == GRETL_OBJ_GPAGE) {
	    display_graph_page();
	} else if (obj->sort == GRETL_OBJ_GRAPH) {
	    open_gui_graph(obj);
	} else if (obj->sort == GRETL_OBJ_PLOT) {
	    open_boxplot(obj);
	}
    } else if (strcmp(item, _("Edit plot commands")) == 0) {
	if (obj->sort == GRETL_OBJ_GRAPH || obj->sort == GRETL_OBJ_PLOT) {
	    GRAPHT *graph = (GRAPHT *) obj->data;

	    remove_png_term_from_plotfile(graph->fname, NULL);
	    view_file(graph->fname, 1, 0, 78, 400, 
		      (obj->sort == GRETL_OBJ_GRAPH)? GR_PLOT : GR_BOX);
	}
    } else if (strcmp(item, _("Delete")) == 0) {
	maybe_delete_session_object(obj);
    } else if (strcmp(item, _("Add to model table")) == 0) {
	if (obj->sort == GRETL_OBJ_EQN) {
	    SESSION_MODEL *mod = (SESSION_MODEL *) obj->data;

	    add_to_model_table(mod->ptr, MODEL_ADD_FROM_MENU, NULL);
	}
    } else if (strcmp(item, _("Clear")) == 0) {
	if (obj->sort == GRETL_OBJ_MODTAB) {
	    clear_model_table(NULL);
	} else if (obj->sort == GRETL_OBJ_GPAGE) {
	    clear_graph_page();
	}
    } else if (strcmp(item, _("Help")) == 0) {
	if (obj->sort == GRETL_OBJ_MODTAB) {
	    context_help(NULL, GINT_TO_POINTER(MODELTAB));
	} else if (obj->sort == GRETL_OBJ_GPAGE) {
	    context_help(NULL, GINT_TO_POINTER(GRAPHPAGE));
	}
    } else if (strcmp(item, _("Options")) == 0) {
	if (obj->sort == GRETL_OBJ_MODTAB) {
	    model_table_dialog();
	} else {
	    dummy_call();
	}
    } else if (strcmp(item, _("Save as TeX...")) == 0) {   
	if (obj->sort == GRETL_OBJ_GPAGE) {
	    graph_page_save_wrapper();
	}
    }
}

static gboolean icon_entered (GtkWidget *icon, GdkEventCrossing *event,
			      gui_obj *gobj)
{
    gtk_widget_set_state(icon, GTK_STATE_SELECTED);
    
    return FALSE;
}

static gboolean icon_left (GtkWidget *icon, GdkEventCrossing *event,
			   gui_obj *gobj)
{
    gtk_widget_set_state(icon, GTK_STATE_NORMAL);
    
    return FALSE;
}

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
    if (info == GRETL_MODEL_POINTER && data != NULL) {
	MODEL **ppmod = (MODEL **) data->data;

	if (ppmod != NULL) {
	    add_to_model_table(*ppmod, MODEL_ADD_BY_DRAG, NULL);
	}
    } else if (info == GRETL_FILENAME && data != NULL) {
	gchar *fname = (gchar *) data->data;

	if (fname != NULL) {
	    graph_page_add_file(fname);
	}
    }
}

static void session_drag_setup (gui_obj *gobj)
{
    GtkWidget *w = GTK_WIDGET(gobj->icon);
    GtkTargetEntry *targ;
    
    if (gobj->sort == GRETL_OBJ_MODTAB) {
	targ = &session_drag_targets[0];
    } else {
	targ = &session_drag_targets[1];
    }

    gtk_drag_dest_set (w,
                       GTK_DEST_DEFAULT_ALL,
                       targ, 1,
                       GDK_ACTION_COPY);

    g_signal_connect (G_OBJECT(w), "drag_data_received",
                      G_CALLBACK(session_data_received),
                      NULL);
}

static void drag_graph (GtkWidget *w, GdkDragContext *context,
			GtkSelectionData *sel, guint info, guint t,
			GRAPHT *graph)
{
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_STRING, 8, 
                           (const guchar *) graph->fname, 
			   strlen(graph->fname));
}

static void graph_drag_connect (GtkWidget *w, GRAPHT *graph)
{
    gtk_drag_source_set(w, GDK_BUTTON1_MASK,
                        &session_drag_targets[1],
                        1, GDK_ACTION_COPY);
    g_signal_connect(G_OBJECT(w), "drag_data_get",
                     G_CALLBACK(drag_graph), graph);
}

static void drag_model (GtkWidget *w, GdkDragContext *context,
			GtkSelectionData *sel, guint info, guint t,
			SESSION_MODEL *mod)
{
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_STRING, 8, 
                           (const guchar *) &mod->ptr, 
			   sizeof mod->ptr);
}

static void model_drag_connect (GtkWidget *w, SESSION_MODEL *mod)
{
    gtk_drag_source_set(w, GDK_BUTTON1_MASK,
                        &session_drag_targets[0],
                        1, GDK_ACTION_COPY);
    g_signal_connect(G_OBJECT(w), "drag_data_get",
                     G_CALLBACK(drag_model), mod);
}

#define WANT_TOOLTIP(t) (t == GRETL_OBJ_EQN || \
	                 t == GRETL_OBJ_GRAPH || \
                         t == GRETL_OBJ_PLOT)

static gui_obj *session_add_icon (gpointer data, int sort, int mode)
{
    gui_obj *gobj;
    gchar *name = NULL;
    GRAPHT *graph = NULL;
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
	graph = (GRAPHT *) data;
	name = g_strdup(graph->name);
	break;
    case GRETL_OBJ_TEXT:
	text = (SESSION_TEXT *) data;
	name = g_strdup(text->name);
	break;
    case GRETL_OBJ_DSET:
	name = g_strdup(_("Data set"));
	break;
    case GRETL_OBJ_INFO:
	name = g_strdup(_("Data info"));
	break;
    case GRETL_OBJ_SCRIPT:
	name = g_strdup(_("Session"));
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
    default:
	break;
    }

    gobj = gui_object_new(name, sort);

    /* full-length object name as tooltip */
    if (strlen(name) > SHOWNAMELEN) {
	gretl_tooltips_add(GTK_WIDGET(gobj->icon), name);
	icon_named = 1;
    }

    /* attach specific data items */
    if (sort == GRETL_OBJ_EQN || sort == GRETL_OBJ_VAR || sort == GRETL_OBJ_SYS) {
	gobj->data = mod;
    } else if (sort == GRETL_OBJ_GRAPH || sort == GRETL_OBJ_PLOT) {
	gobj->data = graph;
    }

    /* set up for drag and drop */
    if (sort == GRETL_OBJ_EQN) {
	model_drag_connect(gobj->icon, gobj->data);
    } else if (sort == GRETL_OBJ_GRAPH) {
	graph_drag_connect(gobj->icon, gobj->data);
    }

    /* second try at adding a tooltip */
    if (WANT_TOOLTIP(sort) && !icon_named) {
	char *str = NULL;
	    
	if (sort == GRETL_OBJ_EQN) {
	    MODEL *pmod = (MODEL *) mod->ptr;

	    str = model_cmd_str(pmod);
	} else if (sort == GRETL_OBJ_GRAPH) {
	    str = graph_str(graph);
	} else if (sort == GRETL_OBJ_PLOT) {
	    str = boxplot_str(graph);
	}
	if (str != NULL) {
	    gretl_tooltips_add(GTK_WIDGET(gobj->icon), str);
	    free(str);
	}
    }

    if      (sort == GRETL_OBJ_TEXT)   gobj->data = text;
    else if (sort == GRETL_OBJ_DSET)   gobj->data = paths.datfile;
    else if (sort == GRETL_OBJ_SCRIPT) gobj->data = cmdfile;
    else if (sort == GRETL_OBJ_MODTAB) gobj->data = NULL;

    if (mode == ICON_ADD_SINGLE) {
	pack_single_icon(gobj);
    } else if (mode == ICON_ADD_BATCH) {
	icon_list = g_list_append(icon_list, gobj);
    }

    return gobj;
}

/* apparatus for rebuilding session models */

static int silent_remember (MODEL **ppmod, DATAINFO *pdinfo)
{
    SESSION_MODEL **models;
    SESSION_MODEL *addmod;
    MODEL *pmod = *ppmod;
    MODEL *tmp; /* will be returned in place of the saved model */
    const char *name = rebuild.model_name[session.nmodels];

#ifdef SESSION_DEBUG
    fprintf(stderr, "silent_remember: session.nmodels = %d\n", session.nmodels);
    fprintf(stderr, "accessing rebuild.model_name[%d]\n", session.nmodels);
#endif

    addmod = session_model_new(pmod, name, GRETL_OBJ_EQN);
    if (addmod == NULL) {
	return 1;
    }

    models = realloc(session.models, (session.nmodels + 1) * sizeof *models);

    if (models == NULL) {
	free_session_model(addmod);
	return 1;
    }

    /* note: augments refcount on pmod */
    gretl_stack_object_as(pmod, GRETL_OBJ_EQN, name);

    session.models = models;
    session.models[session.nmodels] = addmod;
    session.nmodels += 1;

    tmp = gretl_model_new();
    if (tmp == NULL) return 1;

    *ppmod = tmp; /* replaced */

#ifdef SESSION_DEBUG
    fprintf(stderr, "copied '%s' to session.models[%d]\n" 
	    " nmodels = %d\n", name, session.nmodels - 1, 
	    session.nmodels); 
#endif

    return 0;
}

/* The standard behavior here is simply to clear the model for future
   use.  But if we're rebuilding a gretl session then we stack the
   given model pointer and substitute a new blank model for the next
   use.
*/

int clear_or_save_model (MODEL **ppmod, DATAINFO *pdinfo, int rebuilding)
{
    if (rebuilding) {
	static int save;

	if (save) {
	    int i;

#ifdef SESSION_DEBUG
	    fprintf(stderr, "clear_or_save_model: rebuild=%d, model ID=%d,"
		    "rebuild.nmodels=%d\n", rebuilding, (*ppmod)->ID,
		    rebuild.nmodels);
#endif

	    for (i=0; i<rebuild.nmodels; i++) {
		if ((*ppmod)->ID == rebuild.model_ID[i]) {
		    return silent_remember(ppmod, pdinfo);
		}
	    }
	}
	save = 1;
	clear_model(*ppmod);
    } else {
	/* no rebuild, no need to save 
	   FIXME crash here?
	*/
	clear_model(*ppmod);
    }

    return 0;
}

/* delete the file only if it has a generic, non-saved name */

static void maybe_delete_graph_file (const char *fname)
{
    const char *p = strrchr(fname, SLASH);

    if (p != NULL && !strncmp(p, "session.", 8)) {
	remove(fname);
    }
}

void print_saved_object_specs (const char *session_base, FILE *fp)
{
    MODEL *pmod;
    char tmp[MAXLEN];
    int i;

    fprintf(fp, "(* saved objects:\n");

    for (i=0; i<session.nmodels; i++) {
	if (session.models[i]->type == GRETL_OBJ_EQN) {
	    pmod = session.models[i]->ptr;
	    fprintf(fp, "model %d \"%s\"\n", pmod->ID, pmod->name);
	}
    }

    for (i=0; i<session.ngraphs; i++) {
	/* formulate save name for graph */
	sprintf(tmp, "%sGraph_%d", session_base, i + 1);
	/* does the constructed filename (tmp) differ from the
	   current one? */
	if (strcmp(session.graphs[i]->fname, tmp)) {
	    if (copyfile(session.graphs[i]->fname, tmp)) {
		/* copy failed */
		continue;
	    }
	    maybe_delete_graph_file(session.graphs[i]->fname);
	    strcpy(session.graphs[i]->fname, tmp);
	}
	fprintf(fp, "%s %d \"%s\" %s\n", 
		(session.graphs[i]->sort == GRETL_BOXPLOT)?
		"plot" : "graph",
		session.graphs[i]->ID, 
		session.graphs[i]->name, 
		session.graphs[i]->fname);
    }

    fprintf(fp, "*)\n");
}

int print_session_notes (const char *fname)
{
    int err = 0;

    if (session.notes != NULL && strlen(session.notes)) {
	FILE *fp;

	fp = gretl_fopen(fname, "w");
	if (fp != NULL) {
	    fprintf(fp, "%s", session.notes);
	    fclose(fp);
	}  else {
	    err = 1;
	}
    }

    return err;
}

#ifdef OLD_GTK

static GtkWidget *create_popup_item (GtkWidget *popup, char *str, 
				     GtkSignalFunc callback)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(str);
    gtk_signal_connect(GTK_OBJECT(item), "activate",
		       callback,
		       str);
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);

    return item;
}

#else 

static GtkWidget *create_popup_item (GtkWidget *popup, char *str, 
				     GtkCallback callback)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(str);
    g_signal_connect(G_OBJECT(item), "activate",
		     G_CALLBACK(callback),
		     str);
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);

    return item;
}

#endif

static void session_build_popups (void)
{
    size_t i;

    if (global_popup == NULL) {
	global_popup = gtk_menu_new();
	for (i=0; i<sizeof global_items / sizeof global_items[0]; i++) {
	    create_popup_item(global_popup, 
			      _(global_items[i]), 
			      global_popup_activated);
	}
    }

    if (session_popup == NULL) {
	session_popup = gtk_menu_new();
	for (i=0; i<sizeof session_items / sizeof session_items[0]; i++) {
	    create_popup_item(session_popup, 
			      _(session_items[i]), 
			      session_popup_activated);
	}
    }

    if (model_popup == NULL) {
	model_popup = gtk_menu_new();
	for (i=0; i<sizeof model_items / sizeof model_items[0]; i++) {
		create_popup_item(model_popup, 
				  _(model_items[i]), 
				  object_popup_activated);
	}
    }

    if (model_table_popup == NULL) {
	model_table_popup = gtk_menu_new();
	for (i=0; i<sizeof model_table_items / sizeof model_table_items[0]; i++) {
		create_popup_item(model_table_popup, 
				  _(model_table_items[i]), 
				  object_popup_activated);
	}
    }

    if (var_popup == NULL) {
	var_popup = gtk_menu_new();
	for (i=0; i<sizeof var_items / sizeof var_items[0]; i++) {
		create_popup_item(var_popup, 
				  _(var_items[i]), 
				  object_popup_activated);
	}
    }

    if (graph_popup == NULL) {
	graph_popup = gtk_menu_new();
	for (i=0; i<sizeof graph_items / sizeof graph_items[0]; i++) {
		create_popup_item(graph_popup, 
				  _(graph_items[i]), 
				  object_popup_activated);
	}
    }

    if (graph_page_popup == NULL) {
	graph_page_popup = gtk_menu_new();
	for (i=0; i<sizeof graph_page_items / sizeof graph_page_items[0]; i++) {
		create_popup_item(graph_page_popup, 
				  _(graph_page_items[i]), 
				  object_popup_activated);
	}
    }

    if (boxplot_popup == NULL) {
	boxplot_popup = gtk_menu_new();
	for (i=0; i<sizeof graph_items / sizeof graph_items[0]; i++) {
	    if (strstr(graph_items[i], "GUI")) continue;
		create_popup_item(boxplot_popup, 
				  _(graph_items[i]), 
				  object_popup_activated);
	}
    }

    if (data_popup == NULL) {
	data_popup = gtk_menu_new();
	for (i=0; i<sizeof dataset_items / sizeof dataset_items[0]; i++) {
		create_popup_item(data_popup, 
				  _(dataset_items[i]), 
				  data_popup_activated);	    
	}
    }

    if (info_popup == NULL) {
	info_popup = gtk_menu_new();
	for (i=0; i<sizeof info_items / sizeof info_items[0]; i++) {
		create_popup_item(info_popup, 
				  _(info_items[i]), 
				  info_popup_activated);	    
	}
    }
}

static void iconview_connect_signals (GtkWidget *iconview)
{
    g_signal_connect(G_OBJECT(iconview), "destroy",
		     G_CALLBACK(session_view_free), NULL);
    g_signal_connect(G_OBJECT(iconview), "key_press_event",
		     G_CALLBACK(catch_iconview_key), NULL);
}

void view_session (void)
{
    GtkWidget *hbox, *scroller;
    gchar *title;

    if (iconview != NULL) {
	gdk_window_show(iconview->window);
	gdk_window_raise(iconview->window);
	return;
    }

    session_view_init();

    title = g_strdup_printf("gretl: %s", 
			    (session.name[0])? session.name : 
			    _("current session"));
    
    iconview = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(iconview), title);
    g_free(title);
    gtk_window_set_default_size(GTK_WINDOW(iconview), 400, 300);
    gtk_container_set_border_width(GTK_CONTAINER(iconview), 0);

    iconview_connect_signals(iconview);

    session_build_popups();

    hbox = gtk_hbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(iconview), hbox);
    gtk_container_set_border_width(GTK_CONTAINER(hbox), 5);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_container_set_border_width(GTK_CONTAINER(scroller), 0);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(scroller), "button_press_event",
		     G_CALLBACK(session_icon_click), NULL);
#endif

    gtk_box_pack_start(GTK_BOX(hbox), scroller, TRUE, TRUE, 0); 

    icon_table = gtk_table_new(2, SESSION_VIEW_COLS, FALSE);

    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), 
					  icon_table);

    add_all_icons();

    gtk_widget_show(icon_table);
    gtk_widget_show(scroller);
    gtk_widget_show(hbox);
    gtk_widget_show(iconview);

    gtk_container_foreach(GTK_CONTAINER(scroller), 
			 (GtkCallback) white_bg_style, 
			 NULL);

    GTK_WIDGET_SET_FLAGS(icon_table, GTK_CAN_FOCUS);
    gtk_widget_grab_focus(icon_table);
}

/* apparatus for renaming session objects (gtk2 only) */

#ifndef OLD_GTK

static void size_name_entry (GtkWidget *w, const char *name)
{
    PangoLayout *layout;
    PangoRectangle logrect;
    PangoFontDescription *pfd;
#ifdef USE_GNOME
    GtkSettings *settings;
    gchar *fontname;

    settings = gtk_settings_get_default();
    g_object_get(G_OBJECT(settings), "gtk-font-name", &fontname, NULL);
#endif    

    layout = gtk_entry_get_layout(GTK_ENTRY(w));

#ifdef USE_GNOME
    pfd = pango_font_description_from_string(fontname);
    g_free(fontname);
#else
    pfd = pango_font_description_from_string(get_app_fontname());
#endif

    pango_layout_set_font_description(layout, pfd);
    pango_font_description_free(pfd);

    pango_layout_get_pixel_extents(layout, NULL, &logrect);
    gtk_widget_set_size_request(w, 1.1 * logrect.width, -1); 
}

static gboolean object_name_return (GtkWidget *w,
				    GdkEventKey *key,
				    gui_obj *gobj)
{
    if (!gtk_editable_get_editable(GTK_EDITABLE(gobj->label))) {
	return FALSE;
    }

    if (key->keyval == GDK_Return) {
	const gchar *newname = gtk_entry_get_text(GTK_ENTRY(gobj->label));

	gtk_editable_set_position(GTK_EDITABLE(gobj->label), 0);
	gtk_entry_set_has_frame(GTK_ENTRY(gobj->label), FALSE);
	gtk_editable_set_editable(GTK_EDITABLE(gobj->label), FALSE);
	
	if (newname != NULL && *newname != '\0' &&
	    strcmp(newname, gobj->name)) {
	    rename_session_object(gobj, newname);
	}

	size_name_entry(gobj->label, newname);
	gtk_widget_grab_focus(icon_table);

	return TRUE;
    } 

    return FALSE;
}

static gboolean start_rename_object (GtkWidget *w,
				     GdkEventButton *event,
				     gui_obj *gobj)
{
    if (gtk_editable_get_editable(GTK_EDITABLE(gobj->label))) {
	return FALSE;
    }

    gtk_widget_set_size_request(gobj->label, -1, -1);
    gtk_entry_set_width_chars(GTK_ENTRY(gobj->label), SHOWNAMELEN);
    gtk_entry_set_has_frame(GTK_ENTRY(gobj->label), TRUE); 
    gtk_editable_set_editable(GTK_EDITABLE(gobj->label), TRUE);
    gtk_editable_select_region(GTK_EDITABLE(gobj->label), 0, -1);
    gtk_editable_set_position(GTK_EDITABLE(gobj->label), -1);
    gtk_widget_grab_focus(gobj->label);
    
    return TRUE;
}

#endif /* !OLD_GTK -- object renaming stuff */

static void make_short_label_string (char *targ, const char *src)
{
    if (strlen(src) > SHOWNAMELEN) {
	*targ = '\0';
	strncat(targ, src, SHOWNAMELEN - 3);
	strcat(targ, "...");
    } else {
	strcpy(targ, src);
    }
}

#ifdef OLD_GTK

static void create_gobj_icon (gui_obj *gobj, char **xpm)
{
    GdkPixmap *pixmap;
    GdkBitmap *mask;
    GtkStyle *style;
    GtkWidget *image;
    gchar shortname[SHOWNAMELEN + 1];

    style = gtk_widget_get_style(iconview);
    pixmap = gdk_pixmap_create_from_xpm_d(mdata->w->window,
					  &mask, 
					  &style->bg[GTK_STATE_NORMAL], 
					  xpm);

    gobj->icon = gtk_event_box_new();
    gtk_widget_set_usize(gobj->icon, 36, 36);

    image = gtk_pixmap_new(pixmap, mask);

    gtk_container_add(GTK_CONTAINER(gobj->icon), image);
    gtk_widget_show(image);

    if (gobj->sort == GRETL_OBJ_MODTAB || gobj->sort == GRETL_OBJ_GPAGE) {
	session_drag_setup(gobj);
    }

    make_short_label_string(shortname, gobj->name);
    gobj->label = gtk_label_new(shortname);
    
    gtk_signal_connect(GTK_OBJECT(gobj->icon), "button_press_event",
		       GTK_SIGNAL_FUNC(session_icon_click), gobj);
    gtk_signal_connect(GTK_OBJECT(gobj->icon), "enter_notify_event",
		       GTK_SIGNAL_FUNC(icon_entered), gobj);
    gtk_signal_connect(GTK_OBJECT(gobj->icon), "leave_notify_event",
		       GTK_SIGNAL_FUNC(icon_left), gobj);
}

#else

static void create_gobj_icon (gui_obj *gobj, const char **xpm)
{
    GdkPixbuf *pbuf;
    GtkWidget *image;

    pbuf = gdk_pixbuf_new_from_xpm_data(xpm);

    gobj->icon = gtk_event_box_new();
    gtk_widget_set_size_request(gobj->icon, 36, 36);

    image = gtk_image_new_from_pixbuf(pbuf);
    g_object_unref(G_OBJECT(pbuf));

    gtk_container_add(GTK_CONTAINER(gobj->icon), image);
    gtk_widget_show(image);

    if (gobj->sort == GRETL_OBJ_MODTAB || gobj->sort == GRETL_OBJ_GPAGE) {
	session_drag_setup(gobj);
    }

    if (gobj->sort == GRETL_OBJ_EQN || gobj->sort == GRETL_OBJ_GRAPH ||
	gobj->sort == GRETL_OBJ_VAR || gobj->sort == GRETL_OBJ_PLOT || 
	gobj->sort == GRETL_OBJ_SYS) { 
	gobj->label = gtk_entry_new();
	/* on gtk 2.0.N, the text is/was not going into the selected font */
	gtk_entry_set_text(GTK_ENTRY(gobj->label), gobj->name);
	gtk_editable_set_editable(GTK_EDITABLE(gobj->label), FALSE);
	gtk_entry_set_has_frame(GTK_ENTRY(gobj->label), FALSE);
	gtk_entry_set_max_length(GTK_ENTRY(gobj->label), OBJNAMLEN);
	size_name_entry(gobj->label, gobj->name);
	g_signal_connect(G_OBJECT(gobj->label), "button-press-event",
			 G_CALLBACK(start_rename_object), gobj);
	g_signal_connect(G_OBJECT(gobj->label), "key-press-event",
			 G_CALLBACK(object_name_return), gobj);
    } else {
	gchar str[SHOWNAMELEN + 1];

	make_short_label_string(str, gobj->name);
	gobj->label = gtk_label_new(str);
    }

    g_signal_connect(G_OBJECT(gobj->icon), "button_press_event",
		     G_CALLBACK(session_icon_click), gobj);
    g_signal_connect(G_OBJECT(gobj->icon), "enter_notify_event",
		     G_CALLBACK(icon_entered), gobj);
    g_signal_connect(G_OBJECT(gobj->icon), "leave_notify_event",
		     G_CALLBACK(icon_left), gobj);
}

#endif /* OLD_GTK */

static gui_obj *gui_object_new (gchar *name, int sort)
{
    gui_obj *gobj;
    char **xpm = NULL;

    gobj = mymalloc(sizeof *gobj);
    gobj->name = name; 
    gobj->sort = sort;
    gobj->data = NULL;

    switch (sort) {
    case GRETL_OBJ_EQN:
    case GRETL_OBJ_VAR:
    case GRETL_OBJ_SYS:    xpm = model_xpm;       break;
    case GRETL_OBJ_PLOT:   xpm = boxplot_xpm;     break;
    case GRETL_OBJ_GRAPH:  xpm = gnuplot_xpm;     break;
    case GRETL_OBJ_DSET:   xpm = dot_sc_xpm;      break;
    case GRETL_OBJ_INFO:   xpm = xfm_info_xpm;    break;
    case GRETL_OBJ_SCRIPT: xpm = xfm_make_xpm;    break;
    case GRETL_OBJ_TEXT:
    case GRETL_OBJ_NOTES:  xpm = text_xpm;        break;
    case GRETL_OBJ_CORR:   xpm = rhohat_xpm;      break;
    case GRETL_OBJ_STATS:  xpm = summary_xpm;     break;
    case GRETL_OBJ_MODTAB: xpm = model_table_xpm; break;
    case GRETL_OBJ_GPAGE:  xpm = graph_page_xpm;  break;
    default: break;
    }

#ifdef OLD_GTK
    create_gobj_icon(gobj, xpm);
#else
    create_gobj_icon(gobj, (const char **) xpm);
#endif

    return gobj;
} 

#ifdef OLD_GTK

static void auto_save_gp (windata_t *vwin)
{
    FILE *fp;
    gchar *buf;

    if ((fp = gretl_fopen(vwin->fname, "w")) == NULL) {
	buf = g_strdup_printf(_("Couldn't write to %s"), vwin->fname);
	errbox(buf); 
	g_free(buf);
    } else {
	buf = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
	fprintf(fp, "%s", buf);
	g_free(buf); 
	fclose(fp);
    }
}

#else 

static void auto_save_gp (windata_t *vwin)
{
    FILE *fp;
    gchar *buf;
# ifdef ENABLE_NLS
    gchar *trbuf;
# endif

    buf = textview_get_text(GTK_TEXT_VIEW(vwin->w));
    if (buf == NULL) return;

    if ((fp = gretl_fopen(vwin->fname, "w")) == NULL) {
	g_free(buf);
	buf = g_strdup_printf(_("Couldn't write to %s"), vwin->fname);
	errbox(buf); 
	g_free(buf);
	return;
    }

# ifdef ENABLE_NLS
    trbuf = force_locale_from_utf8(buf);
    if (trbuf != NULL) {
	fputs(trbuf, fp);
	g_free(trbuf);
    } else {
	fputs(buf, fp);
    }
# else
    fputs(buf, fp);
# endif

    g_free(buf); 
    fclose(fp);
}

#endif /* OLD_GTK */

#ifdef G_OS_WIN32

static char *add_pause_to_plotfile (const char *fname)
{
    FILE *fin, *fout;
    char fline[MAXLEN];
    char *tmpfile = NULL;
    int gotpause = 0;

    fin = gretl_fopen(fname, "r");
    if (fin == NULL) return NULL;

    tmpfile = g_strdup_printf("%showtmp.gp", paths.userdir);

    fout = gretl_fopen(tmpfile, "w");
    if (fout == NULL) {
	fclose(fin);
	g_free(tmpfile);
	return NULL;
    }

    while (fgets(fline, MAXLEN - 1, fin)) {
	fputs(fline, fout);
	if (strstr(fline, "pause -1")) gotpause = 1;
    }

    if (!gotpause) {
	fputs("pause -1\n", fout);
    }

    fclose(fin);
    fclose(fout);

    return tmpfile;
}

#endif /* G_OS_WIN32 */

#ifdef OLD_GTK

void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w)
{
    gchar *buf = NULL;
    windata_t *vwin = (windata_t *) data;

    auto_save_gp(vwin);

    buf = g_strdup_printf("gnuplot -persist \"%s\"", vwin->fname);

    if (system(buf)) {
        errbox(_("gnuplot command failed"));
    }

    g_free(buf);
}

#else /* gtk-2.0 */

void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w)
{
    gchar *buf = NULL;
    windata_t *vwin = (windata_t *) data;
    int err = 0;
# ifdef G_OS_WIN32
    gchar *tmpfile;
# endif

    auto_save_gp(vwin);

# ifdef G_OS_WIN32
    tmpfile = add_pause_to_plotfile(vwin->fname);
    if (tmpfile != NULL) {
	buf = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, tmpfile);
	err = (WinExec(buf, SW_SHOWNORMAL) < 32);
	remove(tmpfile); /* is this OK? */
	g_free(tmpfile);
    } else {
	err = 1;
    }
# else
    buf = g_strdup_printf("gnuplot -persist \"%s\"", vwin->fname);
    err = gretl_spawn(buf);
# endif

    if (err) {
	errbox(_("gnuplot command failed"));
    }

    g_free(buf);
}

#endif /* gtk versions fork */

void save_plot_commands_callback (GtkWidget *w, gpointer p)
{
    auto_save_gp((windata_t *) p);
}

static void open_gui_graph (gui_obj *gobj)
{
    GRAPHT *graph = (GRAPHT *) gobj->data;

    display_session_graph_png(graph->fname);
}


