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

#ifndef G_OS_WIN32
# include <gtkextra/gtkextra.h>
#else
# include "gtkextra.h"
# include <windows.h>
#endif

#include "pixmaps/model.xpm"
#include "pixmaps/gnuplot.xpm"
#include "pixmaps/xfm_sc.xpm"
#include "pixmaps/xfm_info.xpm"
#include "pixmaps/xfm_text.xpm"
#include "pixmaps/rhohat.xpm"
#include "pixmaps/summary.xpm"

#define SESSION_DEBUG

extern char *endbit (char *dest, char *src, int addscore);
static void auto_save_gp (gpointer data, guint i, GtkWidget *w);
void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w);
gint yes_no_dialog (char *title, char *msg, int cancel);

/* "session" struct and "errtext" are globals */

static char *model_items[] = {
    "Display",
    "Rename",
    "Delete"
};

static char *graph_items[] = {
    "Display",
    "Edit using GUI",
    "Edit plot commands",
    "Rename",
    "Delete"
};

static char *dataset_items[] = {
    "Edit",
    "Save...",
    "Export as CSV..."
};

static char *session_items[] = {
    "Save",
    "Save As...",
    "Add last graph"
};

GtkItemFactoryEntry gp_edit_items[] = {
    { "/_File", NULL, NULL, 0, "<Branch>" }, 
    { "/File/_Save", NULL, auto_save_gp, 0, NULL },
    { "/File/Save _As...", NULL, file_save, SAVE_GP_CMDS, NULL },
    { "/File/Send to _gnuplot", NULL, gp_to_gnuplot, 0, NULL },
    { "/_Edit", NULL, NULL, 0, "<Branch>" },
    { "/Edit/_Copy selection", NULL, text_copy, COPY_SELECTION, NULL },
    { "/Edit/Copy _all", NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

int session_file_open = 0;

static GtkWidget *iconview;
static GtkWidget *session_popup;
static GtkWidget *model_popup;
static GtkWidget *graph_popup;
static GtkWidget *data_popup;
static GtkWidget *iconlist1;
static GList *gobjects;
static GtkIconListItem *active_icon;
static GtkWidget *addgraph;

extern void view_script_default (void);
extern int read_plotfile (GPT_SPEC *plot, char *fname);
extern void show_spreadsheet (char *dataspec);

/* gtkextra functions */
extern void
gtk_icon_list_set_editable (GtkIconList *iconlist, gboolean editable);
extern void
gtk_icon_list_set_text_space (GtkIconList *iconlist, guint spacing);
/* end gtkextra functions */

static void session_build_popups (void);
static void session_popup_activated (GtkWidget *widget, gpointer data);
static void object_popup_activated (GtkWidget *widget, gpointer data);
static void data_popup_activated (GtkWidget *widget, gpointer data);
static void rename_object (GtkWidget *widget, dialog_t *ddata);
static gint check_object_name (gchar *name, gui_obj *gobj, gint sort);

gui_obj *gui_object_new (GtkIconList *iconlist, gchar *name, int sort);
gui_obj *session_new_model (void);
gui_obj *session_add_object (gpointer data, int sort);
void open_gui_model (gui_obj *gobj);
void open_gui_graph (gui_obj *gobj);
static void session_open_object (GtkWidget *widget, 
				 GtkIconListItem *item, GdkEvent *event,
				 gpointer data);

/* ........................................................... */

void session_close_state (gboolean s)
{
    if (mdata->ifac != NULL) 
	flip(mdata->ifac, "/Session/Close", s);
}

/* ........................................................... */

static int rebuild_init (session_t *rebuild)
{
    rebuild->nmodels = 0;
    rebuild->model_ID = malloc(sizeof(int));
    if (rebuild->model_ID == NULL) return 1;
    rebuild->model_name = malloc(sizeof(char *));
    if (rebuild->model_name == NULL) return 1;
    rebuild->model_name[0] = malloc(64);
    if (rebuild->model_name[0] == NULL) return 1;

    return 0;
}

/* ........................................................... */

static void free_rebuild (session_t *rebuild)
{
    int i;

    if (rebuild->model_ID) free(rebuild->model_ID);
    if (rebuild->model_name) {
	for (i=0; i<rebuild->nmodels; i++) {
	    if (rebuild->model_name[i])
		free(rebuild->model_name[i]);
	}
	free(rebuild->model_name);
    }
}

/* .................................................................. */

char *space_to_score (char *str)
{
    int i, n = strlen(str);

    for (i=0; i<n; i++)
	if (str[i] == ' ') str[i] = '_';
    return str;
}

/* ........................................................... */

void add_last_graph (void)
{
    char grname[12], pltname[MAXLEN];
    int i = session.ngraphs;

    sprintf(grname, "Graph %d", plot_count + 1);
    sprintf(pltname, "%ssession.Graph_%d", paths.userdir, plot_count + 1);
    if (copyfile(paths.plotfile, pltname)) {
	errbox("Failed to copy graph commands file");
	return;
    }

    /* write graph into session struct */
    if (session.ngraphs)
	session.graphs = myrealloc(session.graphs, 
				 (i + 1) * sizeof(GRAPHT *));
    else
	session.graphs = mymalloc(sizeof(GRAPHT *));

    if (session.graphs == NULL) return;

    session.graphs[i] = mymalloc(sizeof(GRAPHT));
    if (session.graphs[i] == NULL) return;

    (session.graphs[i])->name = mymalloc(64);
    (session.graphs[i])->fname = mymalloc(64);
    if ((session.graphs[i])->name == NULL || 
	(session.graphs[i])->fname == NULL) {
	return;
    }
    strcpy((session.graphs[i])->fname, pltname);
    strcpy((session.graphs[i])->name, grname);
    (session.graphs[i])->ID = plot_count++;
    session.ngraphs += 1;
    if (iconview == NULL) view_session();
    else session_add_object(session.graphs[i], 'g');   
}

/* ........................................................... */

void remember_model (GtkWidget *widget, dialog_t *ddata)
     /* called from GUI model window, via edit dialog */
{
    windata_t *mydata = ddata->data;
    MODEL *pmod = (MODEL *) mydata->data;
    char *edttext;
    int i = session.nmodels;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (check_object_name(edttext, NULL, 'm')) return;

    if (pmod->name) free(pmod->name);
    if ((pmod->name = mymalloc(64)) == NULL) return;
    strncpy(pmod->name, edttext, 63);
    pmod->name[63] = '\0';

    /* write model into session struct */
    if (session.nmodels)
	session.models = myrealloc(session.models, 
				   (i + 1) * sizeof(MODEL *));
    else
	session.models = mymalloc(sizeof(MODEL *));
    if (session.models == NULL) return;

    session.nmodels += 1;
    session.models[i] = pmod;
    /* session_list(); */
    if (iconview != NULL)
	session_add_object(session.models[i], 'm');
}

/* ........................................................... */

void quick_remember_model (gpointer data, guint j, GtkWidget *widget)
     /* called directly from model window */
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    int i = session.nmodels;

    if (pmod->name) return;
    if ((pmod->name = mymalloc(64)) == NULL) return;
    sprintf(pmod->name, "Model %d", pmod->ID);

    if (session.nmodels)
	session.models = myrealloc(session.models, 
				   (i + 1) * sizeof(MODEL *));
    else
	session.models = mymalloc(sizeof(MODEL *));
    if (session.models == NULL) return;

    session.nmodels += 1;
    session.models[i] = pmod;
    if (iconview != NULL)
	session_add_object(session.models[i], 'm');  

    /* close model window */
    gtk_widget_destroy(gtk_widget_get_toplevel(GTK_WIDGET(mydata->w)));
}

/* ........................................................... */

void session_init (void)
{
    session.models = NULL;
    session.graphs = NULL;
    session.nmodels = 0;
    session.ngraphs = 0;
    session.name[0] = '\0';
}

/* ........................................................... */

void do_open_session (GtkWidget *w, gpointer data)
{
    dialog_t *d = NULL;
    windata_t *fwin = NULL;    

    if (data != NULL) {    
	if (w == NULL) /* not coming from edit_dialog */
	    fwin = (windata_t *) data;
	else {
	    d = (dialog_t *) data;
	    fwin = (windata_t *) d->data;
	}
    }

    clear_data();
    free_session();
    session_init();

#ifdef SESSION_DEBUG
    fprintf(stderr, "do_open_session: about to check %s\n", scriptfile);
#endif

    if (parse_savefile(scriptfile, &session, &rebuild)) 
	return;
    if (recreate_session(scriptfile, &session, &rebuild)) 
	return;

    mkfilelist(2, scriptfile);

    endbit (session.name, scriptfile, 0);

    /* trash the practice files window that launched the query? */
    if (fwin) gtk_widget_destroy(fwin->w);    

    /* sync gui with session */
    session_file_open = 1;
    session_close_state(TRUE);
    session_state(TRUE);
    view_session();
}

/* ........................................................... */

void close_session (void)
{
    clear_data();
    free_session();
    session_state(FALSE);
    session_close_state(FALSE);
    session_file_open = 0;
    /* FIXME - more to do here (icon window?) */
}

/* ........................................................... */

static void free_graph (GRAPHT *graph)
{
    if (graph != NULL) {
	free(graph->name);
	free(graph->fname);
	free(graph);
    }
}

/* ........................................................... */

static int alloc_graph (GRAPHT **pgraph)
{
    *pgraph = mymalloc(sizeof **pgraph);
    if (*pgraph == NULL) return 1;
    (*pgraph)->name = mymalloc(64);
    (*pgraph)->fname = mymalloc(64);
    if ((*pgraph)->name == NULL || (*pgraph)->fname == NULL) 
	return 1;
    return 0;
}

/* ........................................................... */

void free_session (void)
{
    int i;

    if (session.models) {
	for (i=0; i<session.nmodels; i++) 
	    free_model(session.models[i]);
	free(session.models);
    }
    if (session.graphs) {
	for (i=0; i<session.ngraphs; i++) 
	    free_graph(session.graphs[i]);
	free(session.graphs);
    }
    session.nmodels = 0;
    session.models = NULL;
    session.ngraphs = 0;
    session.graphs = NULL;
}

/* ........................................................... */

void print_session (void)
{
    /* testing */
    int i;

    printf("Session contains %d models:\n", session.nmodels);
    for (i=0; i<session.nmodels; i++) {
	printf("model %d %s\n", (session.models[i])->ID,
	       (session.models[i])->name);
    }
    printf("Session contains %d graphs:\n", session.ngraphs);
    for (i=0; i<session.ngraphs; i++) {
	printf("graph: %s (%s)\n", (session.graphs[i])->name,
	       (session.graphs[i])->fname);
    }
}

/* ........................................................... */

int delete_session_model (GtkWidget *w, gpointer data)
{
    dialog_t *myd;
    gui_obj *myobject;
    MODEL *junk, **ppmod;
    int i, j;

    if (data == NULL) return 0; 

    myd = (dialog_t *) data;
    myobject = (gui_obj *) myd->data;
    junk = (MODEL *) myobject->data;
    ppmod = mymalloc((session.nmodels - 1) * sizeof(MODEL *));
    if (session.nmodels > 1 && ppmod == NULL) {
	return 1;
    }
    j = 0;
    for (i=0; i<session.nmodels; i++) {
	if ((session.models[i])->ID != junk->ID) 
	    ppmod[j++] = session.models[i];
	else 
	    free_model(session.models[i]);
    }
    free(session.models);
    session.models = ppmod;
    session.nmodels -= 1;
    gtk_icon_list_remove(GTK_ICON_LIST(iconlist1), active_icon);
    return 0;
}

/* ........................................................... */

int saved_objects (char *fname)
     /* see if there are saved objects from a previous session */
{
    FILE *fp;
    char line[MAXLEN];
    int saves = 0;

#ifdef SESSION_DEBUG
    fprintf(stderr, "saved_objects: checking %s\n", fname);
#endif    

    fp = fopen(fname, "r");
    if (fp == NULL) return -1;

    /* check for objects */
    while (fgets(line, MAXLEN - 1, fp)) {
	if (strncmp(line, "(* saved objects:", 17) == 0) {
	    saves = 1;
	    break;
	}
    }

    fclose(fp);
    return saves;
}

/* ........................................................... */

int parse_savefile (char *fname, SESSION *psession, session_t *rebuild)
{
    FILE *fp;
    char line[MAXLEN], object[7], *tmp;
    int id, i, j, k, n;

    fp = fopen(fname, "r");
    if (fp == NULL) return 1;

    /* find saved objects */
    k = 0;
    while (fgets(line, MAXLEN - 1, fp)) {
	if (strncmp(line, "(* saved objects:", 17) == 0) {
	    k = 1;
	    break;
	}
    }
    if (!k) { /* no saved objects: just a regular script */
	fclose(fp);
	return 1;
    }

    if (rebuild_init(rebuild)) {
	errbox("Out of memory!");
	return 1;
    }

#ifdef SESSION_DEBUG
    fprintf(stderr, "parse_savefile (%s): got saved objects\n", fname);
#endif

    i = 0; /* models */
    k = 0; /* graphs */
    while (fgets(line, MAXLEN - 1, fp)) {
	if (strncmp(line, "*)", 2) == 0) break;
	if (sscanf(line, "%6s %d", object, &id) != 2) {
	    errbox("Session file is corrupted, ignoring");
	    fclose(fp);
	    return 1;
	}
	if (strcmp(object, "model") == 0) {
	    rebuild->nmodels += 1;
#ifdef SESSION_DEBUG
	    fprintf(stderr, "got a model to rebuild (%d)\n"
		    "rebuild->nmodels now = %d\n", id, rebuild->nmodels);
	    
#endif
	    if (i > 0) {
		rebuild->model_ID = myrealloc(rebuild->model_ID,
					      rebuild->nmodels * sizeof(int));
		rebuild->model_name = 
		    myrealloc(rebuild->model_name,
			      rebuild->nmodels * sizeof(char *));
		rebuild->model_name[i] = mymalloc(64);
		if (rebuild->model_ID == NULL ||
		    rebuild->model_name == NULL ||
		    rebuild->model_name[i] == NULL) {
		    fclose(fp);
		    return 1;
		}
	    }
	    rebuild->model_ID[i] = id;
	    tmp = strchr(line, '"') + 1;
	    strncpy(rebuild->model_name[i], tmp, 63);
	    n = strlen(rebuild->model_name[i]);
	    for (j=n; j>0; j--) {
		if (rebuild->model_name[i][j] == '"') {
		    rebuild->model_name[i][j] = '\0';
		    break;
		}
	    }
	    i++;
	    continue;
	}
	if (strcmp(object, "graph") == 0) {
	    psession->ngraphs += 1;
	    if (k > 0) {
		psession->graphs = myrealloc(psession->graphs,
					     psession->ngraphs * 
					     sizeof(GRAPHT *));
	    } else {
		psession->graphs = mymalloc(sizeof(GRAPHT *));
	    }
	    if (psession->graphs == NULL || 
		alloc_graph(&(psession->graphs[k]))) {
		fclose(fp);
		return 1;
	    }
	    tmp = strchr(line, '"') + 1;
	    strncpy((psession->graphs[k])->name, tmp, 63);
	    n = strlen((psession->graphs[k])->name);
	    for (j=n-1; j>0; j--) {
		if ((psession->graphs[k])->name[j] == '"') {
		    (psession->graphs[k])->name[j] = '\0';
		    break;
		}
	    }
	    n = haschar('"', tmp);
	    strcpy((psession->graphs[k])->fname, tmp + n + 1);
	    top_n_tail((psession->graphs[k])->fname);
	    k++;
	    continue;
	}
	else {
	    errbox("Session file is corrupted, ignoring");
	    fclose(fp);
	    return 1;
	}
    }
    fclose(fp);
    return 0;
}

/* ........................................................... */

int recreate_session (char *fname, SESSION *psession, session_t *rebuild)
     /* called on start-up when a "session" file is loaded */
{
    print_t *prn;

    /* no printed output wanted */
    prn = gretl_print_new(GRETL_PRINT_NULL, NULL);

#ifdef SESSION_DEBUG
    fprintf(stderr, "recreate_session: fname = %s\n", fname);
#endif

    if (execute_script(fname, psession, rebuild, prn, REBUILD_EXEC)) 
	errbox("Error recreating session");
    free_rebuild(rebuild);
    gretl_print_destroy(prn);
    return 0;
}

/* ........................................................... */

static void iconview_off (GtkWidget *w, gpointer data)
{
    gui_obj *gobj;

    iconview = NULL;

    while (gobjects != NULL) {
	gobj = (gui_obj *) gobjects->data;
        if (gobj->name) g_free(gobj->name);
	g_free(gobj);
        gobjects = gobjects->next;
    }
    g_list_free(gobjects); /* FIXME? */
}

/* ........................................................... */

static gint null_sort (gpointer a, gpointer b)
{
    return 0;
}

/* ........................................................... */

void view_session (void)
{
    GtkWidget *hbox1, *scrollw1;
    gint i;
    gchar title[80];

    if (iconview != NULL) {
	gdk_window_raise(iconview->window);
	return;
    }

    sprintf(title, "gretl: %s", 
	    (session.name[0])? session.name : "current session");
    
    iconview = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(iconview), title);
    gtk_widget_set_usize(iconview, 400, 300);
    gtk_container_border_width(GTK_CONTAINER(iconview), 0);
    gtk_signal_connect(GTK_OBJECT(iconview), "destroy",
		       GTK_SIGNAL_FUNC(iconview_off), NULL);
    gtk_signal_connect(GTK_OBJECT(iconview), "key_press_event",
		       GTK_SIGNAL_FUNC(catch_key), 
		       (gpointer) iconview);

    session_build_popups();

    hbox1 = gtk_hbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(iconview), hbox1);
    gtk_container_set_border_width(GTK_CONTAINER(hbox1), 5);
    gtk_widget_show(hbox1);

    scrollw1 = gtk_scrolled_window_new(NULL, NULL);
    gtk_container_border_width(GTK_CONTAINER(scrollw1), 0);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollw1),
				   GTK_POLICY_AUTOMATIC,GTK_POLICY_AUTOMATIC);
    gtk_widget_show(scrollw1);

    gtk_box_pack_start(GTK_BOX(hbox1), scrollw1, TRUE, TRUE, 0); 

    iconlist1 = gtk_icon_list_new(48, 2);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollw1),
					  iconlist1);

    gtk_signal_connect(GTK_OBJECT(iconlist1), "select_icon",
		       (GtkSignalFunc) session_open_object, NULL);

    gtk_icon_list_set_editable(GTK_ICON_LIST(iconlist1), FALSE); 
    gtk_icon_list_set_selection_mode(GTK_ICON_LIST(iconlist1), 
				     GTK_SELECTION_BROWSE);
    gtk_icon_list_set_text_space(GTK_ICON_LIST(iconlist1), 80);
    GTK_ICON_LIST(iconlist1)->compare_func = (GCompareFunc) null_sort;

    active_icon = NULL;
    gobjects = NULL;

    if (data_file_open) {
	session_add_object(NULL, 'i');     /* data info */
	session_add_object(NULL, 'd');     /* data file */
	session_add_object(NULL, 'x');     /* summary stats */
	session_add_object(NULL, 'r');     /* correlation matrix */
    }

    session_add_object(NULL, 's');         /* script file */

#ifdef SESSION_DEBUG
    fprintf(stderr, "view_session: session.nmodels = %d\n", session.nmodels);
#endif

    for (i=0; i<session.nmodels; i++) {
#ifdef SESSION_DEBUG
	fprintf(stderr, "adding session.models[%d] to view\n", i);
#endif
	session_add_object(session.models[i], 'm');
    }
    for (i=0; i<session.ngraphs; i++)
	session_add_object(session.graphs[i], 'g');

    gtk_widget_show(iconlist1);
    gtk_widget_show(iconview);
}

/* ........................................................... */

static gint check_object_name (gchar *name, gui_obj *gobj, gint sort)
{
    gint i, n, samename = 0;
    GList *list;
    gui_obj *other;
    gchar *msg;

    if (name == NULL || strlen(name) == 0) {
	errbox("The object must have a name");
	return 1;
    }

    for (n=0; n<strlen(name); n++){
	if(name[n] < 'a' || name[n] > 'z')
	    if(name[n] < 'A' || name[n] > 'Z')
		if(name[n] < '0' || name[n] > '9')
		    if(name[n] != '.' && name[n] != '_' && name[n] != ' '){
			errbox("The name contains invalid characters");
			return 1;
		    }
    }

    if (!gobjects) {
	if (sort == 'm') {
	    for (i=0; i<session.nmodels; i++) 
		if (strcmp(name, (session.models[i])->name) == 0) {
		    samename = 1;
		    break;
		}
	}
	else if (sort == 'g') {
	    for (i=0; i<session.ngraphs; i++) 
		if (strcmp(name, (session.graphs[i])->name) == 0) {
		    samename = 1;
		    break;
		}
	}		    
    }

    list = gobjects;
    while (list) {
	other = (gui_obj *) list->data;
	if (other != gobj
	    && strcmp(other->name, name) == 0
	    && other->sort == sort) {
	    samename = 1;
	    break;
	}
	list = list->next;
    }

    if (samename) {
	msg = g_strdup_printf("Another %s has the same name",
			      (sort == 'm')? "model" : "graph"); 
	errbox(msg);
	g_free(msg);
    }

    return samename;
}

/* ........................................................... */

static void rename_object (GtkWidget *w, dialog_t *ddata)
{
    gchar *name;
    gui_obj *gobj;

    name = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    gobj = (gui_obj *) ddata->data;

    if (gobj->sort != 'm' && gobj->sort != 'g') return;

    if (check_object_name(name, gobj, gobj->sort)) return;

    /* OK, do the actual renaming */
    strcpy(gobj->name, name);
    gtk_icon_list_set_label(GTK_ICON_LIST(iconlist1), active_icon, name);
    if (gobj->sort == 'm') 
	strcpy(((MODEL *) gobj->data)->name, name);
    if (gobj->sort == 'g') 
	strcpy(((GRAPHT *) gobj->data)->name, name);
    
    return;
}

/* ........................................................... */

static void set_addgraph_mode (void)
{
    GtkWidget *gmenu = 
	gtk_item_factory_get_item(mdata->ifac, "/Session/Add last graph");

    gtk_widget_set_sensitive(addgraph, GTK_WIDGET_IS_SENSITIVE(gmenu));
}

/* ........................................................... */

static void object_popup_show (GtkIconListItem *item, GdkEventButton *event)
{
    gui_obj *gobj = (gui_obj *) gtk_icon_list_get_link(item);

    active_icon = item;
    if (gobj->sort == 'm')
	gtk_menu_popup(GTK_MENU(model_popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
    else if (gobj->sort == 'g')
	gtk_menu_popup(GTK_MENU(graph_popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
    else if (gobj->sort == 'd')
	gtk_menu_popup(GTK_MENU(data_popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
    else if (gobj->sort == 's') {
	set_addgraph_mode();
	gtk_menu_popup(GTK_MENU(session_popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
    }
}

/* ........................................................... */

static
void session_open_object (GtkWidget *widget, 
			  GtkIconListItem *item, GdkEvent *event,
			  gpointer data)
{
    GtkIconList *iconlist;
    gui_obj *gobj;
    GdkModifierType mods;
    iconlist = GTK_ICON_LIST(widget);

    gdk_window_get_pointer(widget->window, NULL, NULL, &mods);
    gobj = (gui_obj *) gtk_icon_list_get_link(item);
    if (!event) return;

    if ((mods & GDK_BUTTON1_MASK) && event->type == GDK_2BUTTON_PRESS) {
	switch (gobj->sort) {
	case 'm':
	    open_gui_model(gobj); break;
	case 'g':
	    open_gui_graph(gobj); break;
	case 'd':
	    show_spreadsheet(NULL); break;
	case 'i':
	    open_info(NULL, 0 , NULL); break;
	case 's':
	    view_script_default(); break;
	case 'r':
	    do_menu_op(NULL, CORR, NULL); break;
	case 'x':
	    do_menu_op(NULL, SUMMARY, NULL); break;
	}
    }

    if (mods & GDK_BUTTON3_MASK && 
	(gobj->sort == 'm' || 
	 gobj->sort == 'g' ||
	 gobj->sort == 'd' ||
	 gobj->sort == 's'))
	object_popup_show(item, (GdkEventButton *) event);
  
}


/* ........................................................... */

static void session_build_popups (void)
{
    GtkWidget *item;
    int i;

    if (session_popup == NULL) {
	gint n = sizeof(session_items) / sizeof(session_items[0]);

	session_popup = gtk_menu_new();
	for (i=0; i<n; i++) {
	    if (i < n-1) {
		item = gtk_menu_item_new_with_label(session_items[i]);
		gtk_signal_connect(GTK_OBJECT(item), "activate",
				   (GtkSignalFunc) session_popup_activated,
				   session_items[i]);
		gtk_widget_show(item);
		gtk_menu_append(GTK_MENU(session_popup),item);
	    } else {
		addgraph = gtk_menu_item_new_with_label(session_items[i]);
		gtk_signal_connect(GTK_OBJECT(addgraph), "activate",
				   (GtkSignalFunc) session_popup_activated,
				   session_items[i]);
		gtk_widget_show(addgraph);
		gtk_menu_append(GTK_MENU(session_popup), addgraph);
	    }
	}
    }

    if (model_popup == NULL) {
	model_popup = gtk_menu_new();
	for (i=0; i<sizeof(model_items)/sizeof(model_items[0]); i++) {
	    item = gtk_menu_item_new_with_label(model_items[i]);
	    gtk_signal_connect(GTK_OBJECT(item), "activate",
			       (GtkSignalFunc) object_popup_activated,
			       model_items[i]);
	    gtk_widget_show(item);
	    gtk_menu_append(GTK_MENU(model_popup),item);
	}
    }

    if (graph_popup == NULL) {
	graph_popup = gtk_menu_new();
	for (i=0; i<sizeof(graph_items)/sizeof(graph_items[0]); i++) {
	    item = gtk_menu_item_new_with_label(graph_items[i]);
	    gtk_signal_connect(GTK_OBJECT(item), "activate",
			       (GtkSignalFunc) object_popup_activated,
			       graph_items[i]);
	    gtk_widget_show(item);
	    gtk_menu_append(GTK_MENU(graph_popup),item);
	}
    }

    if (data_popup == NULL) {
	data_popup = gtk_menu_new();
	for (i=0; i<sizeof(dataset_items)/sizeof(dataset_items[0]); i++) {
	    item = gtk_menu_item_new_with_label(dataset_items[i]);
	    gtk_signal_connect(GTK_OBJECT(item), "activate",
			       (GtkSignalFunc) data_popup_activated,
			       dataset_items[i]);
	    gtk_widget_show(item);
	    gtk_menu_append(GTK_MENU(data_popup),item);
	}
    }
}

/* ........................................................... */

static void session_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, "Save") == 0) {
	if (session_file_open && scriptfile[0]) 
	    save_session(scriptfile);
	else
	    file_selector("Save session", paths.userdir, SAVE_SESSION, NULL);
    }
    else if (strcmp(item, "Save As...") == 0) 
	file_selector("Save session", paths.userdir, SAVE_SESSION, NULL);
    else if (strcmp(item, "Add last graph") == 0)
	add_last_graph();
}

/* ........................................................... */

static void data_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, "Edit") == 0) 
	show_spreadsheet(NULL);
    else if (strcmp(item, "Save...") == 0) 
	file_save(mdata, SAVE_DATA, NULL);
    else if (strcmp(item, "Export as CSV...") == 0) 
	file_save(mdata, EXPORT_CSV, NULL);
}

/* ........................................................... */

static void object_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item;
    GtkIconList *iconlist;
    gui_obj *myobject;

    item = (gchar *) data;
    iconlist = GTK_ICON_LIST(iconlist1);

    myobject = (gui_obj *) gtk_icon_list_get_link(active_icon);

    if (strcmp(item, "Display") == 0) {
	if (myobject->sort == 'm') open_gui_model(myobject);
	else if (myobject->sort == 'g') open_gui_graph(myobject);
    } 
    else if (strcmp(item, "Edit using GUI") == 0) {
	if (myobject->sort == 'g') {
	    GRAPHT *graph = (GRAPHT *) myobject->data;
	    GPT_SPEC *plot = mymalloc(sizeof *plot);

	    if (plot == NULL) return;
	    read_plotfile(plot, graph->fname);
	}
    } 
    else if (strcmp(item, "Edit plot commands") == 0) {
	if (myobject->sort == 'g') {
	    GRAPHT *graph = (GRAPHT *) myobject->data;

	    view_file(graph->fname, 1, 0, 78, 400, "gretl: plot commands",
		      gp_edit_items);
	}
    }   
    else if (strcmp(item, "Delete") == 0) {
	gchar text[100];
	sprintf(text,"Really delete '%s'?", myobject->name);
	if (myobject->sort == 'm') {
	    if (!yes_no_dialog("gretl: delete", text, 0)) 
		delete_session_model(NULL, myobject);
	}
    }
    else if (strcmp(item, "Rename") == 0) {
	edit_dialog("gretl: rename", "Enter new name for object:", 
		    myobject->name, 1,
		    "Apply", rename_object, myobject, 
		    "Cancel", NULL, NULL, 0, 0);   
    }
}

/* ........................................................... */

gui_obj *session_add_object (gpointer data, int sort)
{
    gui_obj *gobj;
    gchar *name = NULL;
    MODEL *pmod = NULL;
    GRAPHT *graph = NULL;

    switch (sort) {
    case 'm':
	pmod = (MODEL *) data;
	name = g_strdup(pmod->name);
	break;
    case 'g':
	graph = (GRAPHT *) data;
	name = g_strdup(graph->name);
	break;
    case 'd':
	name = g_strdup("Data set");
	break;
    case 'i':
	name = g_strdup("Data info");
	break;
    case 's':
	name = g_strdup("Session");
	break;
    case 'r':
	name = g_strdup("Corrmat");
	break;
    case 'x':
	name = g_strdup("Summary");
	break;
    default:
	break;
    }

    gobj = gui_object_new(GTK_ICON_LIST(iconlist1), name, sort);
    if (sort == 'm') gobj->data = pmod;
    else if (sort == 'g') gobj->data = graph;
    else if (sort == 'd') gobj->data = paths.datfile;
    else if (sort == 'i') gobj->data = paths.hdrfile;
    else if (sort == 's') gobj->data = cmdfile;
    else if (sort == 'x') gobj->data = NULL;
    g_free(name);

    return gobj;
}

/* ........................................................... */

void open_gui_model (gui_obj *gobj)
{ 
    print_t *prn;
    MODEL *pmod = (MODEL *) gobj->data;

    if (bufopen(&prn)) return;
    printmodel(pmod, datainfo, prn);
    view_model((void *) prn, pmod, 78, 400, gobj->name);
}

/* ........................................................... */

void open_gui_graph (gui_obj *gobj)
{
    char buf[MAXLEN];
    GRAPHT *graph = (GRAPHT *) gobj->data;

#ifdef G_OS_WIN32
    sprintf(buf, "\"%s\" \"%s\"", paths.gnuplot, graph->fname);
    if (WinExec(buf, SW_SHOWNORMAL) < 32)
	errbox("gnuplot command failed");
#else
    sprintf(buf, "gnuplot -persist \"%s\"", graph->fname);
    if (system(buf))
	errbox("gnuplot command failed");
#endif
}

/* ........................................................... */

gui_obj *gui_object_new (GtkIconList *iconlist, gchar *name, int sort)
{
    gui_obj *gobj;
    char **image = NULL;

    gobj = g_new(gui_obj, 1);
    gobj->iconlist = iconlist;
    gobj->name = g_strdup(name); 
    gobj->window = NULL;
    gobj->is_mapped = FALSE;
    gobj->sort = sort;

    gobjects = g_list_append(gobjects, gobj);

    switch (sort) {
    case 'm': image = model_xpm; break;
    case 'g': image = gnuplot_xpm; break;
    case 'd': image = dot_sc_xpm; break;
    case 'i': image = xfm_info_xpm; break;
    case 's': image = text_xpm; break;
    case 'r': image = rhohat_xpm; break;
    case 'x': image = summary_xpm; break;
    default: break;
    }

    gobj->icon = gtk_icon_list_add_from_data(iconlist, image, 
					     name, gobj);
    return gobj;
} 

/* ........................................................... */

static void auto_save_gp (gpointer data, guint quiet, GtkWidget *w)
{
    FILE *fp;
    char msg[MAXLEN];
    gchar *savestuff;
    windata_t *mydata = (windata_t *) data;

    if ((fp = fopen(mydata->fname, "w")) == NULL) {
	sprintf(msg, "couldn't write to %s", mydata->fname);
	errbox(msg); 
	return;
    }
    savestuff = 
	gtk_editable_get_chars(GTK_EDITABLE(mydata->w), 0, -1);
    fprintf(fp, "%s", savestuff);
    g_free(savestuff); 
    fclose(fp);
    if (!quiet) infobox("plot commands saved");
}

/* ........................................................... */

void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w)
{
    char buf[MAXLEN];
    windata_t *mydata = (windata_t *) data;

    auto_save_gp(data, 1, NULL);

#ifdef G_OS_WIN32
    sprintf(buf, "\"%s\" \"%s\"", paths.gnuplot, mydata->fname);
    if (WinExec(buf, SW_SHOWNORMAL) < 32)
        errbox("gnuplot command failed");
#else
    sprintf(buf, "gnuplot -persist \"%s\"", mydata->fname);
    if (system(buf))
        errbox("gnuplot command failed");
#endif
}










