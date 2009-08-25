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
#include "winstack.h"

/* Below: Keep a record of (most) windows that are open, so they can
   be destroyed en masse when a new data file is opened, to prevent
   weirdness that could arise if (e.g.) a model window that pertains
   to a previously opened data file remains open after the data set
   has been changed.  Script windows are exempt, otherwise they are
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
		    mvm = highest_numbered_var_in_model(pmod, datainfo);
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
#if 0
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

    case STACK_MATCH_FNAME:
	if (wstack != NULL) {
	    const char *ctest = (const char *) ptest;

	    for (i=0; i<n_windows; i++) {
		if (wstack[i] != NULL) {
		    windata_t *vwin = 
			g_object_get_data(G_OBJECT(wstack[i]), "vwin");

		    if (vwin != NULL && *vwin->fname != '\0' &&
			strstr(ctest, vwin->fname) != NULL) {
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
    return winstack(STACK_QUERY, NULL, p, NULL);
}

GtkWidget *match_window_by_data (const gpointer p)
{
    GtkWidget *w = NULL;

    winstack(STACK_QUERY, NULL, p, &w);
    return w;
}

GtkWidget *match_window_by_filename (const char *fname)
{
    GtkWidget *w = NULL;

    winstack(STACK_MATCH_FNAME, NULL, fname, &w);
    return w;
}

int highest_numbered_variable_in_winstack (void)
{
    return winstack(STACK_MAXVAR, NULL, NULL, NULL);
}

void winstack_add (GtkWidget *w)
{
    winstack(STACK_ADD, w, NULL, NULL);
}

void winstack_remove (GtkWidget *w)
{
    winstack(STACK_REMOVE, w, NULL, NULL);
}

static void windata_init (windata_t *vwin, int role, gpointer data)
{
    vwin->main = NULL;
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

windata_t *gretl_viewer_new (int role, const gchar *title, 
			     gpointer data, int record)
{
    windata_t *vwin = mymalloc(sizeof *vwin);

    if (vwin == NULL) {
	return NULL;
    }

    windata_init(vwin, role, data);

    if (role == MAINWIN) {
	return vwin;
    }

    vwin->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(vwin->main), title);
    g_signal_connect(G_OBJECT(vwin->main), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    if (record) {
	g_object_set_data(G_OBJECT(vwin->main), "object", data);
	g_object_set_data(G_OBJECT(vwin->main), "role", 
			  GINT_TO_POINTER(vwin->role));
	winstack_add(vwin->main);
    }

    return vwin;
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

    if (record) {
	g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
	g_object_set_data(G_OBJECT(vwin->main), "role", 
			  GINT_TO_POINTER(vwin->role));
	winstack_add(vwin->main);
	g_signal_connect(G_OBJECT(vwin->main), "destroy",
			 G_CALLBACK(free_windata), vwin);
    }	

    return vwin;
}
