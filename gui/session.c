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
#include "selector.h"
#include <sys/stat.h>
#include <unistd.h>

#ifndef G_OS_WIN32
# include <gtkextra/gtkextra.h>
#else
# include "gtkextra.h"
# include <windows.h>
#endif

#include "../pixmaps/model.xpm"
#include "../pixmaps/boxplot.xpm"
#include "../pixmaps/gnuplot.xpm"
#include "../pixmaps/xfm_sc.xpm"
#include "../pixmaps/xfm_info.xpm"
#include "../pixmaps/xfm_text.xpm"
#include "../pixmaps/xfm_make.xpm"
#include "../pixmaps/rhohat.xpm"
#include "../pixmaps/summary.xpm"

/* #define SESSION_DEBUG */

/* from gui_utils.c */
extern void winstack_init (void);
extern void winstack_destroy (void);

static void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w);
static void auto_save_gp (gpointer data, guint i, GtkWidget *w);

/* "session" struct and "errtext" are globals */

static char *model_items[] = {
    N_("Display"),
    N_("Delete")
};

static char *graph_items[] = {
    N_("Display"),
    N_("Edit using GUI"),
    N_("Edit plot commands"),
    N_("Delete")
};

static char *dataset_items[] = {
    N_("Edit"),
    N_("Save..."),
    N_("Export as CSV...")
};

static char *info_items[] = {
    N_("View"),
    N_("Edit"),
};

static char *session_items[] = {
    N_("Save"),
#ifdef GNUPLOT_PNG
    N_("Save As...")
#else
    N_("Save As..."),
    N_("Add last graph")
#endif
};

GtkItemFactoryEntry gp_edit_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" }, 
    { N_("/File/_Save"), NULL, auto_save_gp, 0, NULL },
    { N_("/File/Save _As..."), NULL, file_save, SAVE_GP_CMDS, NULL },
    { N_("/File/Send to _gnuplot"), NULL, gp_to_gnuplot, 0, NULL },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

static int session_file_open;

static GtkWidget *iconview;
static GtkWidget *session_popup;
static GtkWidget *model_popup;
static GtkWidget *graph_popup;
static GtkWidget *data_popup;
static GtkWidget *info_popup;
static GtkWidget *slist;
static GList *gobjects;
static GtkIconListItem *active_icon;
static GtkWidget *addgraph;

extern int read_plotfile (GPT_SPEC *plot, char *fname);
extern void show_spreadsheet (char *dataspec);

/* private functions */
static void session_build_popups (void);
static void session_popup_activated (GtkWidget *widget, gpointer data);
static void object_popup_activated (GtkWidget *widget, gpointer data);
static void data_popup_activated (GtkWidget *widget, gpointer data);
static void info_popup_activated (GtkWidget *widget, gpointer data);
static gui_obj *gui_object_new (GtkIconList *iconlist, gchar *name, int sort);
static gui_obj *session_add_object (gpointer data, int sort);
static void open_gui_model (gui_obj *gobj);
static void open_gui_graph (gui_obj *gobj);
static void open_boxplot (gui_obj *gobj);
static void session_open_object (GtkWidget *widget, 
				 GtkIconListItem *item, GdkEvent *event,
				 gpointer data);

/* ........................................................... */

static int rebuild_init (SESSIONBUILD *rebuild)
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

static void free_rebuild (SESSIONBUILD *rebuild)
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

/* .................................................................. */

static void edit_session_notes (void)
{
    edit_buffer(&session.notes, 80, 400, _("gretl: session notes"),
		EDIT_NOTES);
}

/* .................................................................. */

#ifdef GNUPLOT_PNG

static int filter_plotfile (const char *src, const char *dest)
{
    FILE *fs, *fd;
    char fline[MAXLEN];

    fs = fopen(src, "r");
    if (!fs) {
	errbox("Couldn't read graph file");
	return 1;
    }

    fd = fopen(dest, "w");
    if (!fd) {
	errbox("Couldn't write graph file");
	fclose(fs);
	return 1;
    }

    while (fgets(fline, MAXLEN-1, fs)) {
	if (strncmp(fline, "set term", 8) && 
	    strncmp(fline, "set output", 10))
	    fputs(fline, fd);
    }

    fclose(fs);
    fclose(fd);
    
    return 0;
}

#endif /* GNUPLOT_PNG */

/* ........................................................... */

void add_last_graph (gpointer data, guint code, GtkWidget *w)
{
    char grname[12], pltname[MAXLEN], savedir[MAXLEN];
    int i = session.ngraphs;
    static int boxplot_count;
#ifdef GNUPLOT_PNG
    GPT_SPEC *plot = (GPT_SPEC *) data;
#endif

    get_default_dir(savedir);

    if (code == 0) { /* gnuplot graph */
	sprintf(pltname, "%ssession.Graph_%d", savedir, plot_count + 1);
	sprintf(grname, "%s %d", _("Graph"), plot_count + 1);
#ifdef GNUPLOT_PNG
	if (filter_plotfile(plot->fname, pltname)) return;
#else
	if (copyfile(paths.plotfile, pltname)) {
	    errbox(_("No graph found"));
	    return;
	} 
	remove(paths.plotfile);
#endif
    } else { /* gretl boxplot */
	sprintf(pltname, "%ssession.Plot_%d", savedir, boxplot_count + 1);
	sprintf(grname, "%s %d", _("Boxplot"), boxplot_count + 1);
	boxplot_count++;
	if (copyfile("boxdump.tmp", pltname)) {
	    errbox(_("Failed to copy boxplot file"));
	    return;
	}
	remove("boxdump.tmp");
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

    (session.graphs[i])->sort = code;

    strcpy((session.graphs[i])->fname, pltname);
    strcpy((session.graphs[i])->name, grname);
    (session.graphs[i])->ID = plot_count++;
    session.ngraphs += 1;
    
    session_changed(1);

    if (iconview == NULL)
	infobox(_("Graph saved"));
    else
	session_add_object(session.graphs[i], (code == 0)? 'g' : 'b');
}

/* ........................................................... */

void remember_model (gpointer data, guint close, GtkWidget *widget)
     /* called directly from model window */
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    int i = session.nmodels;
    char buf[24];

    for (i=0; i<session.nmodels; i++) {
	if (session.models[i] == pmod) {
	    infobox(_("Model is already saved"));
	    return;
	}
    }

    pmod->name = g_strdup_printf("%s %d", _("Model"), pmod->ID);

    if (session.nmodels) {
	session.models = myrealloc(session.models, 
				   (i + 1) * sizeof(MODEL *));
    } else {
	session.models = mymalloc(sizeof(MODEL *));
    }

    if (session.models == NULL) return;

    session.nmodels += 1;
    session.models[i] = pmod;

    if (iconview != NULL) {
	session_add_object(session.models[i], 'm'); 
    }

    sprintf(buf, _("%s saved"), pmod->name);
    infobox(buf);

    session_changed(1);

    /* close model window */
    if (close) {
	gtk_widget_destroy(gtk_widget_get_toplevel(GTK_WIDGET(mydata->w)));
    }
}

/* ........................................................... */

int session_changed (int set)
{
    static int has_changed;
    int orig;

    orig = has_changed;
    has_changed = set;
    return orig;
}

/* ........................................................... */

void session_init (void)
{
    session.models = NULL;
    session.graphs = NULL;
    session.notes = NULL;
    session.nmodels = 0;
    session.ngraphs = 0;
    session.name[0] = '\0';
    session_changed(0);
    winstack_init();
}

/* ........................................................... */

void do_open_session (GtkWidget *w, gpointer data)
{
    dialog_t *d = NULL;
    windata_t *fwin = NULL;
    FILE *fp;

    if (data != NULL) {    
	if (w == NULL) /* not coming from edit_dialog */
	    fwin = (windata_t *) data;
	else {
	    d = (dialog_t *) data;
	    fwin = (windata_t *) d->data;
	}
    }

    fp = fopen(tryscript, "r");
    if (fp != NULL) {
	fclose(fp);
	strcpy(scriptfile, tryscript);
    } else {
	char errbuf[MAXLEN];

	sprintf(errbuf, _("Couldn't open %s\n"), tryscript);
	errbox(errbuf);
	delete_from_filelist(2, tryscript);
	delete_from_filelist(3, tryscript);
	return;
    }

    clear_data(1);
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

    endbit(session.name, scriptfile, 0);

    /* pick up session notes, if any */
    if (1) {
	char notesfile[MAXLEN];
	struct stat buf;

	switch_ext(notesfile, scriptfile, "Notes");
	if (stat(notesfile, &buf) == 0 && (fp = fopen(notesfile, "r"))) { 
	    char notesline[MAXLEN];

	    /* read into buffer */
	    session.notes = mymalloc(buf.st_size + 1);
	    if (session.notes == NULL) {
		fclose(fp);
		return;
	    }
	    session.notes[0] = '\0';
	    while (fgets(notesline, MAXLEN-1, fp)) {
		strcat(session.notes, notesline);
	    } 
	    fclose(fp);
	}
    }

    /* trash the practice files window that launched the query? */
    if (fwin) gtk_widget_destroy(fwin->w);    

    /* sync gui with session */
    session_file_open = 1;
    session_state(TRUE);
    view_session();
}

/* ........................................................... */

void close_session (void)
{
    clear_data(1); /* this was (0): why?? */
    free_session();
    session_state(FALSE);
    session_file_open = 0;
    if (iconview != NULL) 
	gtk_widget_destroy(iconview);
    session_changed(0);
    winstack_destroy();
    clear_selector();
}

/* ........................................................... */

void verify_clear_data (void)
{
    if (!expert) {
        int button = yes_no_dialog ("gretl",                      
				    _("Clearing the data set will end\n"
				    "your current session.  Continue?"), 0);
        if (button != YES_BUTTON) 
            return;
    }
    close_session();
}

/* ........................................................... */

void free_session (void)
{
    int i;

    if (session.models) {
	for (i=0; i<session.nmodels; i++) 
	    free_model(session.models[i]);
	free(session.models);
	session.models = NULL;
    }
    if (session.graphs) {
	free(session.graphs);
	session.graphs = NULL;
    }
    if (session.notes) {
	free(session.notes);
	session.notes = NULL;
    }
    session.nmodels = 0;
    session.ngraphs = 0;
    session.name[0] = '\0';
}

/* ........................................................... */

#ifdef SESSION_DEBUG
void print_session (void)
{
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
#endif

/* ........................................................... */

static int delete_session_object (gui_obj *obj)
{
    int i, j;

    if (obj == NULL) return 0; 

    if (obj->sort == 'm') { /* it's a model */
	MODEL **ppmod, *junk = (MODEL *) obj->data;

	/* special case: only one model currently */
	if (session.nmodels == 1) {
	    free_model(session.models[0]);
	} else {
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
	}
	session.nmodels -= 1;
    }
    else if (obj->sort == 'g') { /* it's a graph */    
	GRAPHT **ppgr, *junk = (GRAPHT *) obj->data;

	/* special case: only one graph currently */
	if (session.ngraphs == 1) {
	    free(session.graphs[0]);
	} else {
	    ppgr = mymalloc((session.ngraphs - 1) * sizeof(GRAPHT *));
	    if (session.ngraphs > 1 && ppgr == NULL) {
		return 1;
	    }
	    j = 0;
	    for (i=0; i<session.ngraphs; i++) {
		if ((session.graphs[i])->ID != junk->ID) 
		    ppgr[j++] = session.graphs[i];
		else 
		    free(session.graphs[i]);
	    }
	    free(session.graphs);
	    session.graphs = ppgr;
	}
	session.ngraphs -= 1;
    }

    session_changed(1);

    gtk_icon_list_remove(GTK_ICON_LIST(slist), active_icon);
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

int parse_savefile (char *fname, SESSION *psession, SESSIONBUILD *rebuild)
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
	errbox(_("Out of memory!"));
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
	    errbox(_("Session file is corrupted, ignoring"));
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
		rebuild->model_name[i] = mymalloc(24);
		if (rebuild->model_ID == NULL ||
		    rebuild->model_name == NULL ||
		    rebuild->model_name[i] == NULL) {
		    fclose(fp);
		    return 1;
		}
	    }
	    rebuild->model_ID[i] = id;
	    tmp = strchr(line, '"') + 1;
	    strncpy(rebuild->model_name[i], tmp, 23);
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
	if (!strcmp(object, "graph") || !strcmp(object, "plot")) {
	    psession->ngraphs += 1;
	    if (k > 0) {
		psession->graphs = myrealloc(psession->graphs,
					     psession->ngraphs * 
					     sizeof(GRAPHT *));
	    } else {
		psession->graphs = mymalloc(sizeof(GRAPHT *));
	    }
	    if (psession->graphs == NULL) {
		fclose(fp);
		return 1;
	    }
	    psession->graphs[k] = mymalloc(sizeof(GRAPHT));
	    if (psession->graphs[k] == NULL) {
		fclose(fp);
		return 1;
	    }
	    tmp = strchr(line, '"') + 1;
	    strncpy((psession->graphs[k])->name, tmp, 23);
	    (psession->graphs[k])->name[23] = '\0';
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
#ifdef SESSION_DEBUG
	    fprintf(stderr, "got graph: '%s'\n", (psession->graphs[k])->fname);
#endif
	    (psession->graphs[k])->ID = plot_count++;
	    k++;
	    continue;
	} else {
	    errbox(_("Session file is corrupted, ignoring"));
	    fclose(fp);
	    return 1;
	}
    }
#ifdef SESSION_DEBUG
    fprintf(stderr, "psession->ngraphs = %d\n", psession->ngraphs);
#endif
    fclose(fp);
    return 0;
}

/* ........................................................... */

int recreate_session (char *fname, SESSION *psession, SESSIONBUILD *rebuild)
     /* called on start-up when a "session" file is loaded */
{
    PRN *prn;
    extern int replay; /* lib.c */

    /* no printed output wanted */
    prn = gretl_print_new(GRETL_PRINT_NULL, NULL);

#ifdef SESSION_DEBUG
    fprintf(stderr, "recreate_session: fname = %s\n", fname);
    print_session();
#endif

    if (execute_script(fname, NULL, psession, rebuild, prn, REBUILD_EXEC)) 
	errbox(_("Error recreating session"));
#ifdef SESSION_DEBUG
    fprintf(stderr, "recreate_session: after execute_script()\n");
    print_session();
#endif
    free_rebuild(rebuild);
    gretl_print_destroy(prn);
    replay = 1; /* no fresh commands have been entered yet */
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
	gdk_window_show(iconview->window);
	gdk_window_raise(iconview->window);
	return;
    }

    sprintf(title, "gretl: %s", 
	    (session.name[0])? session.name : _("current session"));
    
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

    slist = gtk_icon_list_new(48, 2);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollw1),
					  slist);

    gtk_signal_connect(GTK_OBJECT(slist), "select_icon",
		       (GtkSignalFunc) session_open_object, NULL);

    gtk_icon_list_set_editable(GTK_ICON_LIST(slist), FALSE); 
    gtk_icon_list_set_selection_mode(GTK_ICON_LIST(slist), 
				     GTK_SELECTION_BROWSE);
    gtk_icon_list_set_text_space(GTK_ICON_LIST(slist), 80);
    GTK_ICON_LIST(slist)->compare_func = (GCompareFunc) null_sort;

    active_icon = NULL;
    gobjects = NULL;

    if (data_status) {
	session_add_object(NULL, 'i');     /* data info */
	session_add_object(NULL, 'd');     /* data file */
	session_add_object(NULL, 'n');     /* session notes */
	session_add_object(NULL, 'x');     /* summary stats */
	session_add_object(NULL, 'r');     /* correlation matrix */
    }

    session_add_object(NULL, 's');         /* script file */

#ifdef SESSION_DEBUG
    fprintf(stderr, "view_session: session.nmodels = %d\n", session.nmodels);
    fprintf(stderr, "view_session: session.ngraphs = %d\n", session.ngraphs);
#endif

    for (i=0; i<session.nmodels; i++) {
	
#ifdef SESSION_DEBUG
	fprintf(stderr, "adding session.models[%d] to view\n", i);
#endif
	session_add_object(session.models[i], 'm');
    }
    for (i=0; i<session.ngraphs; i++) {
#ifdef SESSION_DEBUG
	fprintf(stderr, "adding session.graphs[%d] to view\n", i);
#endif
	/* distinguish gnuplot graphs from gretl boxplots */
	session_add_object(session.graphs[i], 
			   ((session.graphs[i])->sort == 1)? 'b' : 'g');
    }

    gtk_widget_show(slist);
    gtk_widget_show(iconview);
}

/* ........................................................... */

static void set_addgraph_mode (void)
{
    GtkWidget *gmenu = 
	gtk_item_factory_get_item(mdata->ifac, _("/Session/Add last graph"));

    if (gmenu == NULL || addgraph == NULL) return;

    gtk_widget_set_sensitive(addgraph, GTK_WIDGET_IS_SENSITIVE(gmenu));
}

/* ........................................................... */

static void object_popup_show (GtkIconListItem *item, GdkEventButton *event)
{
    gui_obj *gobj = (gui_obj *) gtk_icon_list_get_link(item);
    GtkWidget *w = NULL;

    active_icon = item;

    switch (gobj->sort) {
    case 'm': 
	w = model_popup; 
	break;
    case 'g': 
	w = graph_popup; 
	break;
    case 'd': 
	w = data_popup; 
	break;
    case 'i': 
	w = info_popup; 
	break;
    case 's': 
	set_addgraph_mode(); 
	w = session_popup; 
	break;
    default: 
	break;
    }
	
    gtk_menu_popup(GTK_MENU(w), NULL, NULL, NULL, NULL,
		   event->button, event->time);
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
	case 'b':
	    open_boxplot(gobj); break;
	case 'g':
	    open_gui_graph(gobj); break;
	case 'd':
	    show_spreadsheet(NULL); break;
	case 'i':
	    open_info(NULL, 0, NULL); break;
	case 's':
	    view_script_default(); break;
	case 'n':
	    edit_session_notes(); break;
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
         gobj->sort == 'i' ||
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
		item = gtk_menu_item_new_with_label(_(session_items[i]));
		gtk_signal_connect(GTK_OBJECT(item), "activate",
				   (GtkSignalFunc) session_popup_activated,
				   _(session_items[i]));
		gtk_widget_show(item);
		gtk_menu_append(GTK_MENU(session_popup),item);
	    } else {
		addgraph = gtk_menu_item_new_with_label(_(session_items[i]));
		gtk_signal_connect(GTK_OBJECT(addgraph), "activate",
				   (GtkSignalFunc) session_popup_activated,
				   _(session_items[i]));
		gtk_widget_show(addgraph);
		gtk_menu_append(GTK_MENU(session_popup), addgraph);
	    }
	}
    }

    if (model_popup == NULL) {
	model_popup = gtk_menu_new();
	for (i=0; i<sizeof(model_items)/sizeof(model_items[0]); i++) {
	    item = gtk_menu_item_new_with_label(_(model_items[i]));
	    gtk_signal_connect(GTK_OBJECT(item), "activate",
			       (GtkSignalFunc) object_popup_activated,
			       _(model_items[i]));
	    gtk_widget_show(item);
	    gtk_menu_append(GTK_MENU(model_popup),item);
	}
    }

    if (graph_popup == NULL) {
	graph_popup = gtk_menu_new();
	for (i=0; i<sizeof(graph_items)/sizeof(graph_items[0]); i++) {
	    item = gtk_menu_item_new_with_label(_(graph_items[i]));
	    gtk_signal_connect(GTK_OBJECT(item), "activate",
			       (GtkSignalFunc) object_popup_activated,
			       _(graph_items[i]));
	    gtk_widget_show(item);
	    gtk_menu_append(GTK_MENU(graph_popup),item);
	}
    }

    if (data_popup == NULL) {
	data_popup = gtk_menu_new();
	for (i=0; i<sizeof(dataset_items)/sizeof(dataset_items[0]); i++) {
	    item = gtk_menu_item_new_with_label(_(dataset_items[i]));
	    gtk_signal_connect(GTK_OBJECT(item), "activate",
			       (GtkSignalFunc) data_popup_activated,
			       _(dataset_items[i]));
	    gtk_widget_show(item);
	    gtk_menu_append(GTK_MENU(data_popup),item);
	}
    }

    if (info_popup == NULL) {
	info_popup = gtk_menu_new();
	for (i=0; i<sizeof(info_items)/sizeof(info_items[0]); i++) {
	    item = gtk_menu_item_new_with_label(_(info_items[i]));
	    gtk_signal_connect(GTK_OBJECT(item), "activate",
			       (GtkSignalFunc) info_popup_activated,
			       _(info_items[i]));
	    gtk_widget_show(item);
	    gtk_menu_append(GTK_MENU(info_popup),item);
	}
    }
}

/* ........................................................... */

void save_session_callback (GtkWidget *w, guint i, gpointer data)
{
    if (i == 0 && session_file_open && scriptfile[0]) {
	save_session(scriptfile);
	session_changed(0);
    } else {
	file_selector(_("Save session"), SAVE_SESSION, NULL);
    }
}

/* ........................................................... */

static void session_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("Save")) == 0) 
	save_session_callback(NULL, 0, NULL);
    else if (strcmp(item, _("Save As...")) == 0) 
	save_session_callback(NULL, 1, NULL);
    else if (strcmp(item, _("Add last graph")) == 0)
	add_last_graph(NULL, 0, NULL);
}

/* ........................................................... */

static void info_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("View")) == 0) 
	open_info(NULL, 0, NULL);
    else if (strcmp(item, _("Edit")) == 0) 
	edit_header(NULL, 0, NULL);
}

/* ........................................................... */

static void data_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (strcmp(item, _("Edit")) == 0) 
	show_spreadsheet(NULL);
    else if (strcmp(item, _("Save...")) == 0) 
	file_save(mdata, SAVE_DATA, NULL);
    else if (strcmp(item, _("Export as CSV...")) == 0) 
	file_save(mdata, EXPORT_CSV, NULL);
}

/* ........................................................... */

static void object_popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item;
    GtkIconList *iconlist;
    gui_obj *myobject;

    item = (gchar *) data;
    iconlist = GTK_ICON_LIST(slist);

    myobject = (gui_obj *) gtk_icon_list_get_link(active_icon);

    if (strcmp(item, _("Display")) == 0) {
	if (myobject->sort == 'm') open_gui_model(myobject);
	else if (myobject->sort == 'g') open_gui_graph(myobject);
    } 
    else if (strcmp(item, _("Edit using GUI")) == 0) {
	if (myobject->sort == 'g') {
	    GRAPHT *graph = (GRAPHT *) myobject->data;
	    GPT_SPEC *plot = mymalloc(sizeof *plot);

	    if (plot == NULL) return;
	    read_plotfile(plot, graph->fname);
	}
    } 
    else if (strcmp(item, _("Edit plot commands")) == 0) {
	if (myobject->sort == 'g') {
	    GRAPHT *graph = (GRAPHT *) myobject->data;

	    view_file(graph->fname, 1, 0, 78, 400, GR_PLOT, gp_edit_items);
	}
    }   
    else if (strcmp(item, _("Delete")) == 0) {
	gchar msg[64];

	sprintf(msg, _("Really delete %s?"), myobject->name);
	if (!yes_no_dialog(_("gretl: delete"), msg, 0)) {
	    delete_session_object(myobject);
	}
    }
}

/* ........................................................... */

static void store_list (int *list, char *buf)
{
    int i;
    char numstr[5];

    for (i=1; i<=list[0]; i++) {
        sprintf(numstr, "%d ", list[i]);
        strcat(buf, numstr);
    }
}

/* ........................................................... */

static char *model_cmd_str (MODEL *pmod)
{
    char *str;
    
    str = malloc(MAXLEN);
    if (str == NULL) return NULL;

    sprintf(str, "%s ", commands[pmod->ci]);

    if (pmod->ci == AR) {
        store_list(pmod->arinfo->arlist, str);
        strcat(str, "; ");
    }
    store_list(pmod->list, str);    

    return str;
}

/* ........................................................... */

static char *graph_str (GRAPHT *graph)
{
    FILE *fp;
    char *str = NULL;

    fp = fopen(graph->fname, "r");
    if (fp != NULL) {
	char xlabel[24], ylabel[24], line[48];
	int gotxy = 0;

	while (fgets(line, 47, fp) && gotxy < 2) {
	    if (sscanf(line, "set xlabel %23s", xlabel) == 1) 
		gotxy++;
	    else if (sscanf(line, "set ylabel %23s", ylabel) == 1)
		gotxy++;
	}
	if (gotxy == 2 && (str = malloc(64))) {
	    sprintf(str, "%s versus %s", ylabel, xlabel);
	}
	fclose(fp);
    }
    return str;
}

/* ........................................................... */

static char *boxplot_str (GRAPHT *graph)
{
    FILE *fp;
    char *str = NULL;

    fp = fopen(graph->fname, "r");
    if (fp != NULL) {
	char vname[9], line[48];

	str = malloc (MAXLEN);
	if (str == NULL) return NULL;
	str[0] = '\0';

	while (fgets(line, 47, fp) && strlen(str) < MAXLEN-10) {
	    if (sscanf(line, "%*d varname = %8s", vname) == 1) { 
		strcat(str, vname);
		strcat(str, " ");
	    }
	}
	fclose(fp);
    }
    return str;
}

/* ........................................................... */

static gui_obj *session_add_object (gpointer data, int sort)
{
    gui_obj *gobj;
    gchar *name = NULL;
    MODEL *pmod = NULL;
    GRAPHT *graph = NULL;
    extern GtkTooltips *gretl_tips;

    switch (sort) {
    case 'm':
	pmod = (MODEL *) data;
	name = g_strdup(pmod->name);
	break;
    case 'b':
    case 'g':
	graph = (GRAPHT *) data;
	name = g_strdup(graph->name);
	break;
    case 'd':
	name = g_strdup(_("Data set"));
	break;
    case 'i':
	name = g_strdup(_("Data info"));
	break;
    case 's':
	name = g_strdup(_("Session"));
	break;
    case 'n':
	name = g_strdup(_("Notes"));
	break;
    case 'r':
	name = g_strdup(_("Corrmat"));
	break;
    case 'x':
	name = g_strdup(_("Summary"));
	break;
    default:
	break;
    }

    gobj = gui_object_new(GTK_ICON_LIST(slist), name, sort);
    if (sort == 'm') {
	char *str = model_cmd_str(pmod);

	gobj->data = pmod;
	if (str != NULL) {
	    gtk_tooltips_set_tip(gretl_tips, GTK_WIDGET(gobj->icon->entry), 
				 str, NULL);
	    free(str);
	}
    }
    else if (sort == 'g') {
	char *str = graph_str(graph);

	gobj->data = graph;
	if (str != NULL) {
	    gtk_tooltips_set_tip(gretl_tips, GTK_WIDGET(gobj->icon->entry), 
				 str, NULL);
	    free(str);
	}
    }    
    else if (sort == 'b') {
	char *str = boxplot_str(graph);

	gobj->data = graph;
	if (str != NULL) {
	    gtk_tooltips_set_tip(gretl_tips, GTK_WIDGET(gobj->icon->entry), 
				 str, NULL);
	    free(str);
	}
    }
    else if (sort == 'd') gobj->data = paths.datfile;
    else if (sort == 'i') gobj->data = NULL;
    else if (sort == 's') gobj->data = cmdfile;
    else if (sort == 'n') gobj->data = NULL;
    else if (sort == 'x') gobj->data = NULL;
    g_free(name);

    return gobj;
}

/* ........................................................... */

static void open_gui_model (gui_obj *gobj)
{ 
    PRN *prn;
    MODEL *pmod = (MODEL *) gobj->data;

    if (bufopen(&prn)) return;
    if (printmodel(pmod, datainfo, prn))
	pmod->errcode = E_NAN;
    view_model((void *) prn, pmod, 78, 400, gobj->name);
}

/* ........................................................... */

static void open_gui_graph (gui_obj *gobj)
{
    char buf[MAXLEN];
    GRAPHT *graph = (GRAPHT *) gobj->data;

#ifdef G_OS_WIN32
    sprintf(buf, "\"%s\" \"%s\"", paths.gnuplot, graph->fname);
    if (WinExec(buf, SW_SHOWNORMAL) < 32)
	errbox(_("gnuplot command failed"));
#else
    sprintf(buf, "gnuplot -persist \"%s\"", graph->fname);
    if (system(buf))
	errbox(_("gnuplot command failed"));
#endif
}

/* ........................................................... */

extern int retrieve_boxplot (const char *fname);

static void open_boxplot (gui_obj *gobj)
{
    GRAPHT *graph = (GRAPHT *) gobj->data;

    if (retrieve_boxplot(graph->fname)) 
	errbox(_("Failed to reconstruct boxplot"));
}

/* ........................................................... */

static gui_obj *gui_object_new (GtkIconList *iconlist, gchar *name, int sort)
{
    gui_obj *gobj;
    char **image = NULL;

    gobj = g_new(gui_obj, 1);
    gobj->name = g_strdup(name); 
    gobj->sort = sort;

    gobjects = g_list_append(gobjects, gobj);

    switch (sort) {
    case 'm': image = model_xpm; break;
    case 'b': image = boxplot_xpm; break;
    case 'g': image = gnuplot_xpm; break;
    case 'd': image = dot_sc_xpm; break;
    case 'i': image = xfm_info_xpm; break;
    case 's': image = xfm_make_xpm; break;
    case 'n': image = text_xpm; break;
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
	sprintf(msg, _("Couldn't write to %s"), mydata->fname);
	errbox(msg); 
	return;
    }
    savestuff = 
	gtk_editable_get_chars(GTK_EDITABLE(mydata->w), 0, -1);
    fprintf(fp, "%s", savestuff);
    g_free(savestuff); 
    fclose(fp);
    if (!quiet) infobox(_("plot commands saved"));
}

/* ........................................................... */

static void gp_to_gnuplot (gpointer data, guint i, GtkWidget *w)
{
    char buf[MAXLEN];
    windata_t *mydata = (windata_t *) data;

    auto_save_gp(data, 1, NULL);

#ifdef G_OS_WIN32
    sprintf(buf, "\"%s\" \"%s\"", paths.gnuplot, mydata->fname);
    if (WinExec(buf, SW_SHOWNORMAL) < 32)
        errbox(_("gnuplot command failed"));
#else
    sprintf(buf, "gnuplot -persist \"%s\"", mydata->fname);
    if (system(buf))
        errbox(_("gnuplot command failed"));
#endif
}










