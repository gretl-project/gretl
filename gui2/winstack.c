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

#include "gretl.h"
#include "var.h"
#include "dlgutils.h"
#include "winstack.h"

#define WDEBUG 0

/* Below: Keep a record of (most) windows that are open, so they can
   be destroyed en masse when a new data file is opened, to prevent
   weirdness that could arise if (e.g.) a model window that pertains
   to a previously opened data file remains open after the data set
   has been changed. Script windows are exempt, otherwise they are
   likely to disappear when their "run" control is activated, which we
   don't want.
*/

enum winstack_codes {
    STACK_INIT,
    STACK_ADD,
    STACK_REMOVE,
    STACK_DESTROY,
    STACK_QUERY,
    STACK_MATCH_FNAME,
    STACK_MATCH_FNAME_MOD,
    STACK_MATCH_VWIN,
    STACK_MAXVAR
};

static int max_var_in_stacked_models (GtkWidget **wstack, int nwin)
{
    int i, role, mvm, vmax = 0;

    for (i=0; i<nwin; i++) {
	if (wstack[i] != NULL) {
	    role = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(wstack[i]), "role"));
	    if (role == VIEW_MODEL) {
		const MODEL *pmod;

		pmod = g_object_get_data(G_OBJECT(wstack[i]), "object");
		if (pmod != NULL) {
		    mvm = highest_numbered_var_in_model(pmod, dataset);
		    if (mvm > vmax) {
			vmax = mvm;
		    }
		}
	    } else if (role == VAR || role == VECM) {
		const GRETL_VAR *var;

		var = g_object_get_data(G_OBJECT(wstack[i]), "object");
		if (var != NULL) {
		    mvm = gretl_VAR_get_highest_variable(var);
		    if (mvm > vmax) {
			vmax = mvm;
		    }		    
		}
	    } 
	}
    }    

    return vmax;
}

static int got_filename_match (GtkWidget *w, const char *s,
			       int code)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(w), "vwin");

    if (vwin != NULL && *vwin->fname != '\0') {
	if (code == STACK_MATCH_FNAME) {
	    return !strcmp(s, vwin->fname);
	} else {
	    /* vwin->fname may have the suffix removed */
	    return !strncmp(s, vwin->fname, strlen(vwin->fname));
	}
    } else {
	return 0;
    }
}

static int 
winstack (int code, GtkWidget *w, gconstpointer ptest, GtkWidget **pw)
{
    static int n_windows;
    static GtkWidget **wstack;
    int i, ret = 0;

    switch (code) {

    case STACK_DESTROY:	
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] != NULL) {
#if WDEBUG
		fprintf(stderr, "winstack: destroying widget at %p\n", 
			(void *) wstack[i]);
#endif
		gtk_widget_destroy(wstack[i]);
	    }
	}
	free(wstack);
	/* fall-through intended */

    case STACK_INIT:
	wstack = NULL;
	n_windows = 0;
	break;

    case STACK_ADD:
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] == NULL) {
		wstack[i] = w;
		break;
	    }
	}
	if (i == n_windows) {
	    GtkWidget **newstack;

	    newstack = myrealloc(wstack, (n_windows + 1) * sizeof *wstack);
	    if (newstack != NULL) { 
		wstack = newstack;
		wstack[n_windows] = w;
		n_windows++;
	    }
	}
	break;

    case STACK_REMOVE:
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] == w) {
		wstack[i] = NULL;
		ret = 1;
		break;
	    }
	}
	break;

    case STACK_QUERY:
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] != NULL) {
		gpointer p = g_object_get_data(G_OBJECT(wstack[i]), "object");

		if (p == ptest) {
		    if (pw != NULL) {
			*pw = wstack[i];
		    }
		    ret = 1;
		    break;
		}
	    }
	}
	break;

    case STACK_MATCH_VWIN:
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] == ptest) {
		ret = 1;
		break;
	    }
	}
	break;

    case STACK_MATCH_FNAME:
    case STACK_MATCH_FNAME_MOD:
	if (wstack != NULL) {
	    const char *ctest = (const char *) ptest;

	    for (i=0; i<n_windows; i++) {
		if (wstack[i] != NULL) {
		    if (got_filename_match(wstack[i], ctest, code)) {
			if (pw != NULL) {
			    *pw = wstack[i];
			}
			ret = 1;
			break;
		    }
		}
	    }
	}
	break;

    case STACK_MAXVAR:
	ret = max_var_in_stacked_models(wstack, n_windows);
	break;	

    default:
	break;
    }

    return ret;
}

void winstack_init (void)
{
    winstack(STACK_INIT, NULL, NULL, NULL);
}
    
void winstack_destroy (void)
{
    winstack(STACK_DESTROY, NULL, NULL, NULL);
}

int winstack_match_data (const gpointer p)
{
    GtkWidget *w = NULL;
    int ret = winstack(STACK_QUERY, NULL, p, &w);
    
    if (w != NULL) {
	gtk_window_present(GTK_WINDOW(w));
    }

    return ret;
}

int vwin_on_stack (const windata_t *vwin)
{
    return winstack(STACK_MATCH_VWIN, NULL, vwin->main, NULL);
}

GtkWidget *match_window_by_data (const gpointer p)
{
    GtkWidget *w = NULL;

    winstack(STACK_QUERY, NULL, p, &w);
    return w;
}

void maybe_close_window_for_data (const gpointer p,
				  GretlObjType otype)
{
    GtkWidget *w = NULL;

    winstack(STACK_QUERY, NULL, p, &w);
    if (w != NULL) {
	if (otype == GRETL_OBJ_BUNDLE) {
	    windata_t *vwin = g_object_get_data(G_OBJECT(w),
						"vwin");
	    if (vwin != NULL) {
		vwin->data = NULL;
	    }
	} else if (otype == GRETL_OBJ_MATRIX) {
	    g_object_set_data(G_OBJECT(w), "object", NULL);
	}
	gtk_widget_destroy(w);
    }
}

GtkWidget *match_window_by_filename (const char *fname)
{
    GtkWidget *w = NULL;

    winstack(STACK_MATCH_FNAME, NULL, fname, &w);
    return w;
}

GtkWidget *match_db_window_by_filename (const char *fname)
{
    GtkWidget *w = NULL;

    winstack(STACK_MATCH_FNAME_MOD, NULL, fname, &w);
    return w;
}

int highest_numbered_variable_in_winstack (void)
{
    return winstack(STACK_MAXVAR, NULL, NULL, NULL);
}

void winstack_add (GtkWidget *w)
{
#if WDEBUG
    fprintf(stderr, "winstack add: %p (%s)\n", (void *) w,
	    gtk_window_get_title(GTK_WINDOW(w)));
#endif
    winstack(STACK_ADD, w, NULL, NULL);
}

void winstack_remove (GtkWidget *w)
{
#if WDEBUG
    fprintf(stderr, "winstack remove: %p (%s)\n", (void *) w,
	    gtk_window_get_title(GTK_WINDOW(w)));    
#endif
    winstack(STACK_REMOVE, w, NULL, NULL);
}

static void windata_init (windata_t *vwin, int role, gpointer data)
{
    vwin->main = NULL;
    vwin->topmain = NULL;
    vwin->vbox = NULL;
    vwin->text = NULL;
    vwin->listbox = NULL;
    vwin->mbar = NULL;
    vwin->finder = NULL;
    vwin->status = NULL;
    vwin->popup = NULL;
    vwin->ui = NULL;
    vwin->gretl_parent = NULL;
    vwin->gretl_children = NULL;
    vwin->data = data;
    vwin->active_var = 0;
    vwin->role = role;
    vwin->n_model_tests = 0;
    vwin->n_gretl_children = 0;
    vwin->flags = 0;
    vwin->fname[0] = '\0';
    vwin->sbuf = NULL;
}

static void viewer_box_config (windata_t *vwin)
{
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 0);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);
    
#ifndef G_OS_WIN32
    set_wm_icon(vwin_toplevel(vwin));
#endif
}

windata_t *
gretl_viewer_new_with_parent (windata_t *parent, int role, 
			      const gchar *title, 
			      gpointer data, int record)
{
    windata_t *vwin = mymalloc(sizeof *vwin);

    if (vwin == NULL) {
	return NULL;
    }

    windata_init(vwin, role, data);

    if (role == MAINWIN) {
	/* gets special treatment in gretl.c */
	return vwin;
    }

    vwin->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(vwin->main), title);
    g_signal_connect(G_OBJECT(vwin->main), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);

    if (record) {
	g_object_set_data(G_OBJECT(vwin->main), "object", data);
	g_object_set_data(G_OBJECT(vwin->main), "role", 
			  GINT_TO_POINTER(vwin->role));
	winstack_add(vwin->main);
    } 

    if (parent != NULL) {
	vwin_add_child(parent, vwin);
    }

    viewer_box_config(vwin);
    add_window_list_item(vwin->main, role);

    return vwin;
}

static void free_tabwin (tabwin_t *twin)
{
    return;
}

static void destroy_viewer_tab (GtkWidget *w, windata_t *vwin)
{
    gtk_widget_destroy(vwin->main);
}

static GtkWidget *
make_viewer_tab (windata_t *vwin, const gchar *filename)
{
    gchar *title = NULL;
    const gchar *p;
    GtkWidget *tab;
    GtkWidget *label;
    GtkWidget *img;
    GtkWidget *button;

    tab = gtk_hbox_new(FALSE, 0);

    if (strstr(filename, "script_tmp") != NULL) {
	title = g_strdup("untitled");
    } else if ((p = strrchr(filename, SLASH)) != NULL) {
	title = g_strdup(p + 1);
    } else {
	title = g_strdup(filename);
    }

    label = gtk_label_new(title);
    gtk_container_add(GTK_CONTAINER(tab), label);
    g_free(title);

    img = gtk_image_new_from_stock(GTK_STOCK_CLOSE, 
				   GTK_ICON_SIZE_MENU);
    button = gtk_button_new();
    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    gtk_container_add(GTK_CONTAINER(button), img);
    g_signal_connect(button, "clicked", G_CALLBACK(destroy_viewer_tab), 
		     vwin);

    gtk_container_add(GTK_CONTAINER(tab), button);
    gtk_widget_show_all(tab);
    
    return tab;
}

windata_t *gretl_tabbed_viewer_new (int role, const gchar *title, 
				    const gchar *filename, 
				    gpointer data, int record)
{
    tabwin_t *twin = mymalloc(sizeof *twin);
    windata_t *vwin;
    GtkWidget *hbox;
    GtkWidget *tab;

    if (twin == NULL) {
	return NULL;
    }

    vwin = mymalloc(sizeof *vwin);
    if (vwin == NULL) {
	free(twin);
	return NULL;
    }

    windata_init(vwin, role, data);
    vwin->main = gtk_hbox_new(FALSE, 0);
    g_signal_connect(G_OBJECT(vwin->main), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    twin->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(twin->main), title);
    g_signal_connect(G_OBJECT(twin->main), "destroy", 
		     G_CALLBACK(free_tabwin), vwin);
    g_object_set_data(G_OBJECT(twin->main), "twin", twin);

    /* vertically oriented container */
    twin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(twin->vbox), 0);
    gtk_box_set_spacing(GTK_BOX(twin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(twin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(twin->main), twin->vbox);

    /* hbox to hold menu bar */
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(twin->vbox), hbox, FALSE, FALSE, 0);
    g_object_set_data(G_OBJECT(twin->main), "hbox", hbox);

    twin->tabs = gtk_notebook_new();
    gtk_container_add(GTK_CONTAINER(twin->vbox), twin->tabs);

    tab = make_viewer_tab(vwin, filename);
    gtk_notebook_append_page(GTK_NOTEBOOK(twin->tabs), vwin->main, tab);
    vwin->topmain = twin->main;

    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);

#if 0 /* FIXME */
    if (record) {
	g_object_set_data(G_OBJECT(vwin->main), "object", data);
	g_object_set_data(G_OBJECT(vwin->main), "role", 
			  GINT_TO_POINTER(vwin->role));
	winstack_add(vwin->main);
    } 
#endif

    viewer_box_config(vwin);
    add_window_list_item(twin->main, role);

    return vwin;
}

windata_t *gretl_viewer_new (int role, const gchar *title, 
			     gpointer data, int record)
{
    return gretl_viewer_new_with_parent(NULL, role, title,
					data, record);
}

GtkWidget *vwin_toplevel (windata_t *vwin)
{
    return vwin->topmain != NULL ? vwin->topmain : vwin->main;
}

void vwin_pack_toolbar (windata_t *vwin)
{
    GtkWidget *hbox;

    if (vwin->topmain != NULL) {
	hbox = g_object_get_data(G_OBJECT(vwin->topmain), "hbox");
	gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    } else {    
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
	gtk_widget_show_all(hbox);
    }
}

static gint catch_winlist_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    GdkModifierType mods = widget_get_pointer_mask(w);

    if ((mods & GDK_MOD1_MASK) && key->keyval == GDK_w) {
	window_list_popup(w, NULL, vwin->main);
	return TRUE;
    }

    return FALSE;
}

windata_t *gretl_browser_new (int role, const gchar *title, int record)
{
    windata_t *vwin = mymalloc(sizeof *vwin);

    if (vwin == NULL) {
	return NULL;
    }

    windata_init(vwin, role, NULL);

    vwin->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(vwin->main), title);

    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_winlist_key), vwin);

    if (record) {
	g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
	g_object_set_data(G_OBJECT(vwin->main), "role", 
			  GINT_TO_POINTER(vwin->role));
	winstack_add(vwin->main);
	g_signal_connect(G_OBJECT(vwin->main), "destroy",
			 G_CALLBACK(free_windata), vwin);
    } 

    add_window_list_item(vwin->main, role);

    return vwin;
}
