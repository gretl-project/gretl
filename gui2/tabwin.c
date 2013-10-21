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
#include "winstack.h"
#include "textbuf.h"
#include "tabwin.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#define TDEBUG 0

struct tabwin_t_ {
    int role;             /* what's tabwin doing? */ 
    GtkWidget *main;      /* top-level GTK window */
    GtkWidget *mbox;      /* horizontal box to hold menu bar */
    GtkWidget *mbar;      /* menu bar */
    GtkWidget *tabs;      /* notebook for tabs */
    GtkWidget *dialog;    /* associated dialog */
    GtkWidget *dlg_owner; /* the tab that "owns" dialog */
};

GtkTargetEntry tabwin_drag_targets[] = {
    { "text/uri-list",  0, GRETL_FILENAME },
};

/* We support a tabbed editor for gretl scripts and also
   a tabbed viewer for models */

static tabwin_t *tabedit;
static tabwin_t *tabmod;

static void tabwin_destroy (GtkWidget *w, tabwin_t *tabwin)
{
    if (tabwin == tabedit) {
	tabedit = NULL;
    } else if (tabwin == tabmod) {
	tabmod = NULL;
    }

    free(tabwin);
}

static tabwin_t *vwin_get_tabwin (windata_t *vwin)
{
    return g_object_get_data(G_OBJECT(vwin->topmain), "tabwin");
}

static gboolean maybe_block_tabedit_quit (tabwin_t *tabwin,
					  GtkWindow *topwin)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    int np = gtk_notebook_get_n_pages(notebook);
    gboolean ret = FALSE;

    if (tabwin->role != EDIT_SCRIPT) {
	return FALSE;
    }

    if (np > 1) {
	gchar *msg = g_strdup_printf(_("Editing %d scripts: really quit?"), np);
	gint resp;

	if (topwin != NULL) {
	    gtk_window_present(topwin);
	}

	resp = yes_no_dialog(_("gretl: script editor"), msg, 0);
	if (resp != GRETL_YES) {
	    ret = TRUE;
	}
	g_free(msg);
    } else if (np == 1) {
	GtkWidget *page = gtk_notebook_get_nth_page(notebook, 0);

	if (page != NULL) {
	    windata_t *vwin = g_object_get_data(G_OBJECT(page), "vwin");
    
	    if (vwin_content_changed(vwin)) {
		if (topwin != NULL) {
		    gtk_window_present(topwin);
		}
		ret = query_save_text(NULL, NULL, vwin);
	    }
	}
    }

    return ret;
}

/* called on program exit */

gboolean tabwin_exit_check (GtkWidget *w)
{
    tabwin_t *tabwin = g_object_get_data(G_OBJECT(w), "tabwin");

    return maybe_block_tabedit_quit(tabwin, GTK_WINDOW(w));
}

/* called on delete-event */

static gboolean tabedit_quit_check (GtkWidget *w, GdkEvent *event, 
				    tabwin_t *tabwin)
{
    return maybe_block_tabedit_quit(tabwin, NULL);
}

/* called by top-level toolbar closer button */

void maybe_destroy_tabwin (windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);

    if (!maybe_block_tabedit_quit(tabwin, NULL)) {
	gtk_widget_destroy(vwin->topmain);
    }
}

/* activate or de-activate a tab's closer button */

static void viewer_tab_show_closer (GtkNotebook *notebook,
				    GtkWidget *tab, 
				    gboolean show)
{
    GtkWidget *lbl = gtk_notebook_get_tab_label(notebook, tab);    
    GtkWidget *button = g_object_get_data(G_OBJECT(lbl), "closer");

    if (button != NULL) {
	if (show) {
	    gtk_widget_show(button);
	} else {
	    gtk_widget_hide(button);
	}
    }
}

static void tabwin_remove_toolbar (tabwin_t *tabwin)
{
#if TDEBUG
    fprintf(stderr, "*** removing toolbar at %p\n", (void *) tabwin->mbar);
#endif
    gtk_container_remove(GTK_CONTAINER(tabwin->mbox), 
			 tabwin->mbar);
    tabwin->mbar = NULL;
}

static void tabwin_insert_toolbar (tabwin_t *tabwin, windata_t *vwin)
{
#if TDEBUG
    fprintf(stderr, "*** inserting toolbar at %p\n", (void *) vwin->mbar);
#endif
    gtk_box_pack_start(GTK_BOX(tabwin->mbox), vwin->mbar, 
		       TRUE, TRUE, 0);
    gtk_widget_show_all(vwin->mbar);
    tabwin->mbar = vwin->mbar;
}

static void page_removed_callback (GtkNotebook *notebook,
				   GtkWidget *child,
				   gint pgnum,
				   gpointer data)
{
    int np = gtk_notebook_get_n_pages(notebook);

#if TDEBUG
    fprintf(stderr, "*** page_removed_callback: child=%p\n", (void *) child);
#endif

    if (np < 5) {
	gtk_notebook_popup_disable(notebook);
    }    

    if (np == 1) {
	/* only one tab left after removal: this page should
	   not display its own closer button, nor should it
	   be detachable
	*/
	GtkWidget *tab = gtk_notebook_get_nth_page(notebook, 0);

	gtk_notebook_set_tab_detachable(notebook, tab, FALSE);
	viewer_tab_show_closer(notebook, tab, FALSE);
    }
}

static void page_added_callback (GtkNotebook *notebook,
				 GtkWidget *child,
				 gint pgnum,
				 gpointer data)
{
    int i, np = gtk_notebook_get_n_pages(notebook);
    GtkWidget *tab;

    if (np >= 5) {
	gtk_notebook_popup_enable(notebook);
    }

    if (np > 1) {
	for (i=0; i<np; i++) {
	    tab = gtk_notebook_get_nth_page(notebook, i);
	    gtk_notebook_set_tab_detachable(notebook, tab, TRUE);
	    viewer_tab_show_closer(notebook, tab, TRUE);
	}
    } else {
	tab = gtk_notebook_get_nth_page(notebook, 0);
	viewer_tab_show_closer(notebook, tab, FALSE);
    }
}

/* callback for tab-specific close button: this should be
   invoked only if there's at least one tab left after
   trashing the selected one
*/

static void tabwin_tab_close (GtkWidget *w, windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    gint pg = gtk_notebook_page_num(notebook, vwin->main);

#if TDEBUG
    fprintf(stderr, "*** tabwin_tab_close: vwin = %p\n", (void *) vwin);
#endif

    if (vwin->main == tabwin->dlg_owner) {
	/* this tab has a dialog open: don't close it */
	gtk_window_present(GTK_WINDOW(tabwin->dialog));
	return;
    }

    /* note: vwin->mbar is packed under tabwin, so it will not
       get destroyed automatically when the page is removed
    */
    if (tabwin->mbar != NULL && tabwin->mbar == vwin->mbar) {
	tabwin_remove_toolbar(tabwin);
    }

    /* relinquish the extra reference */
#if TDEBUG
    fprintf(stderr, " unrefing toolbar at %p\n", (void *) vwin->mbar);
#endif
    g_object_unref(vwin->mbar);

    gtk_notebook_remove_page(notebook, pg);
}

void tabwin_tab_destroy (windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);

#if TDEBUG
    fprintf(stderr, "*** tabwin_tab_destroy: vwin = %p\n", (void *) vwin);
#endif
    
    if (gtk_notebook_get_n_pages(notebook) > 1) {
	gint pg = gtk_notebook_page_num(notebook, vwin->main);

	if (tabwin->mbar != NULL && tabwin->mbar == vwin->mbar) {
	    tabwin_remove_toolbar(tabwin);
	}
	g_object_unref(vwin->mbar);
	gtk_notebook_remove_page(notebook, pg);
    } else {
	gtk_widget_destroy(vwin->topmain);
    }
}

/* on switching the current page, put the new page's
   toolbar into place in tabwin (and remove the old
   one, if present)
*/

static gboolean switch_page_callback (GtkNotebook *tabs,
				      gpointer arg1,
				      gint pgnum,
				      tabwin_t *tabwin)
{
    GtkWidget *tab = gtk_notebook_get_nth_page(tabs, pgnum);
    windata_t *vwin = NULL;

    if (tab != NULL) {
	vwin = g_object_get_data(G_OBJECT(tab), "vwin");
    }

#if TDEBUG
    fprintf(stderr, "*** switch_page_callback: tab=%p, vwin=%p\n",
	    (void *) tab, (void *) vwin);
#endif

    if (vwin != NULL) {
	if (tabwin->mbar != NULL && tabwin->mbar != vwin->mbar) {
	    /* there's an "old" toolbar in place */
	    tabwin_remove_toolbar(tabwin);
	}
	if (vwin->mbar != NULL && vwin->mbar != tabwin->mbar) {
	    /* a "new" toolbar should be shown */
	    tabwin_insert_toolbar(tabwin, vwin);
	}
    }    

    return FALSE;
}

/* callback for the "create-window" signal */

static GtkNotebook *detach_tab_callback (GtkNotebook *book,
					 GtkWidget *page,
					 gint x, gint y,
					 gpointer data)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(page), "vwin");

    if (vwin != NULL) {
	undock_tabbed_viewer(NULL, vwin);
    }

    return NULL;
}

/* avoid excessive padding in a tab's close button */

static void no_button_padding (GtkWidget *w)
{
    static int style_done;

    gtk_widget_set_name(w, "closer");

    if (!style_done) {
	gtk_rc_parse_string("style \"closer-style\"\n{\n"
			    "  GtkWidget::focus-padding = 0\n"
			    "  GtkWidget::focus-line-width = 0\n"
			    "  xthickness = 0\n"
			    "  ythickness = 0\n"
			    "}\n"
			    "widget \"*.closer\" style \"closer-style\"");
	style_done = 1;
    }
}

/* put a tab-specific close button next to the tab's label */

static void viewer_tab_add_closer (GtkWidget *tab, windata_t *vwin)
{
    GtkWidget *img, *button;
	
    img = gtk_image_new_from_stock(GTK_STOCK_CLOSE, 
				   GTK_ICON_SIZE_MENU);
    button = gtk_button_new();
    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    gtk_container_set_border_width(GTK_CONTAINER(button), 0);
    no_button_padding(button);
    gtk_container_add(GTK_CONTAINER(button), img);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(tabwin_tab_close), 
		     vwin);
    gtk_container_add(GTK_CONTAINER(tab), button);
    g_object_set_data(G_OBJECT(tab), "closer", button);
}

/* try to ensure unique dummy title strings for unsaved
   new scripts */

static gchar *untitled_title (tabwin_t *tabwin)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    GtkWidget *page;
    windata_t *vwin;
    int i, np, idx = 0;

    np = gtk_notebook_get_n_pages(notebook);

    for (i=0; i<np; i++) {
	page = gtk_notebook_get_nth_page(notebook, i);
	vwin = g_object_get_data(G_OBJECT(page), "vwin");
	if (strstr(vwin->fname, "script_tmp") != NULL) {
	    idx++;
	}
    }

    if (idx > 0) {
	return g_strdup_printf("%s(%d)", _("untitled"), idx);
    } else {
	return g_strdup(_("untitled"));
    }
}

/* create and add tab with filename and closer button */

static GtkWidget *make_viewer_tab (tabwin_t *tabwin, 
				   windata_t *vwin, 
				   const gchar *info)
{
    GtkNotebook *notebook;
    gchar *title = NULL;
    GtkWidget *tab;
    GtkWidget *label;
    GtkWidget *mlabel;

    tab = gtk_hbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(tab), 0);

    if (tabwin->role == EDIT_SCRIPT) {
	if (strstr(info, "script_tmp") != NULL) {
	    title = untitled_title(tabwin);
	} else {
	    title = title_from_filename(info, FALSE);
	}
    } else if (info != NULL) {
	title = g_strdup(info);
    } else {
	title = g_strdup("unknown");
    }

    label = gtk_label_new(title);
    mlabel = gtk_label_new(title);
    gtk_container_add(GTK_CONTAINER(tab), label);
    g_object_set_data(G_OBJECT(tab), "label", label);
    g_object_set_data(G_OBJECT(tab), "mlabel", mlabel);
    g_free(title);

    viewer_tab_add_closer(tab, vwin);
    gtk_widget_set_size_request(tab, -1, 18);
    gtk_widget_show_all(tab);

    notebook = GTK_NOTEBOOK(tabwin->tabs);
    gtk_notebook_append_page_menu(notebook, vwin->main, 
				  tab, mlabel);

    return tab;
}

static gint catch_tabwin_key (GtkWidget *w, GdkEventKey *key, 
			      tabwin_t *tabwin)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    gint pg = gtk_notebook_get_current_page(notebook);
    GtkWidget *tab = gtk_notebook_get_nth_page(notebook, pg);
    windata_t *vwin = g_object_get_data(G_OBJECT(tab), "vwin");

    return catch_viewer_key(w, key, vwin);
}

static void  
tabwin_handle_drag  (GtkWidget *widget,
		     GdkDragContext *context,
		     gint x,
		     gint y,
		     GtkSelectionData *data,
		     guint info,
		     guint time,
		     gpointer p)
{
    const guchar *seldata = NULL;
    gchar *dfname;
    char tmp[MAXLEN];
    int pos, skip = 5;

    if (data != NULL) {
	seldata = gtk_selection_data_get_data(data);
    }

    if (info != GRETL_FILENAME) {
	return;
    }

    /* ignore the wrong sort of data */
    if (data == NULL || (dfname = (gchar *) seldata) == NULL || 
	strlen(dfname) <= 5 || strncmp(dfname, "file:", 5)) {
	return;
    }

    if (strncmp(dfname, "file://", 7) == 0) skip = 7;
#ifdef G_OS_WIN32
    if (strncmp(dfname, "file:///", 8) == 0) skip = 8;
#endif

    /* there may be multiple files: we ignore all but the first */
    *tmp = 0;
    if ((pos = gretl_charpos('\r', dfname)) > 0 || 
	(pos = gretl_charpos('\n', dfname) > 0)) {
	strncat(tmp, dfname + skip, pos - skip);
    } else {
	strcat(tmp, dfname + skip);
    }

    /* handle spaces and such */
    unescape_url(tmp);

#ifdef G_OS_WIN32
    filename_to_win32(tryfile, tmp);
#else
    strcpy(tryfile, tmp);
#endif

    if (has_suffix(tryfile, ".inp")) {
	do_open_script(EDIT_SCRIPT);
    }
}

/* build a tabbed viewer/editor */

static tabwin_t *make_tabbed_viewer (int role)
{
    tabwin_t *tabwin;
    GtkWidget *vbox;

    tabwin = mymalloc(sizeof *tabwin);

    if (tabwin == NULL) {
	return NULL;
    }

    tabwin->role = role;
    tabwin->dialog = NULL;
    tabwin->dlg_owner = NULL;

    /* top-level window */
    tabwin->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    if (role == EDIT_SCRIPT) {
	gtk_window_set_title(GTK_WINDOW(tabwin->main), 
			     _("gretl: script editor"));
 	g_signal_connect(G_OBJECT(tabwin->main), "delete-event",
			 G_CALLBACK(tabedit_quit_check), tabwin);
    } else {
	gtk_window_set_title(GTK_WINDOW(tabwin->main), _("gretl: models"));
    }	
    g_signal_connect(G_OBJECT(tabwin->main), "destroy", 
		     G_CALLBACK(tabwin_destroy), tabwin);
    g_object_set_data(G_OBJECT(tabwin->main), "tabwin", tabwin);

    /* vertically oriented container */
    vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vbox), 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 0);
    gtk_container_add(GTK_CONTAINER(tabwin->main), vbox);

    /* box to hold menu bar */
    tabwin->mbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), tabwin->mbox, FALSE, FALSE, 0);
    tabwin->mbar = NULL;

    /* notebook with its signal handlers */
    tabwin->tabs = gtk_notebook_new();
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(tabwin->tabs), TRUE);
    g_signal_connect(G_OBJECT(tabwin->tabs), "switch-page",
		     G_CALLBACK(switch_page_callback), tabwin);
    g_signal_connect(G_OBJECT(tabwin->tabs), "create-window",
		     G_CALLBACK(detach_tab_callback), tabwin);
    g_signal_connect(G_OBJECT(tabwin->tabs), "page-added",
		     G_CALLBACK(page_added_callback), tabwin);
    g_signal_connect(G_OBJECT(tabwin->tabs), "page-removed",
		     G_CALLBACK(page_removed_callback), tabwin);
    gtk_container_add(GTK_CONTAINER(vbox), tabwin->tabs);

#ifndef G_OS_WIN32
    set_wm_icon(tabwin->main);
#endif

    return tabwin;
}

static tabwin_t *get_tabwin_for_role (int role, int *starting)
{
    tabwin_t *tabwin = NULL;

    if (role == EDIT_SCRIPT) {
	if (tabedit != NULL) {
	    tabwin = tabedit;
	} else {
	    *starting = 1;
	    tabedit = tabwin = make_tabbed_viewer(role);
	}
    } else if (role == VIEW_MODEL) {
	if (tabmod != NULL) {
	    tabwin = tabmod;
	} else {
	    *starting = 1;
	    tabmod = tabwin = make_tabbed_viewer(role);
	}
    }	

    return tabwin;
}

/* Create a viewer/editor tab, as an alternative to a stand-alone
   window. We build the tabbed top-level if need be, otherwise we
   stick a new tab into the existing window.
*/

windata_t *viewer_tab_new (int role, const char *info,
			   gpointer data)
{
    tabwin_t *tabwin;
    windata_t *vwin;
    int starting = 0;
    gulong handler_id;

    tabwin = get_tabwin_for_role(role, &starting);
    if (tabwin == NULL) {
	return NULL;
    }

    vwin = vwin_new(role, data);
    if (vwin == NULL) {
	return NULL;
    }

    vwin->flags = VWIN_TABBED;
    vwin->main = gtk_hbox_new(FALSE, 0);
    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
    handler_id = g_signal_connect(G_OBJECT(vwin->main), "destroy", 
				  G_CALLBACK(free_windata), vwin);
    g_object_set_data(G_OBJECT(vwin->main), "destroy-id",
		      GUINT_TO_POINTER(handler_id));

#if TDEBUG
    fprintf(stderr, "*** viewer_tab_new: vwin=%p, main hbox=%p\n", 
	    (void *) vwin, (void *) vwin->main);
#endif

    make_viewer_tab(tabwin, vwin, info);
    vwin->topmain = tabwin->main;

    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 1);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    if (starting) {
	window_list_add(tabwin->main, role);
	g_signal_connect(G_OBJECT(tabwin->main), "key-press-event", 
			 G_CALLBACK(catch_tabwin_key), tabwin);
    } 

    return vwin;
}

/* called when a new editor tab is added: if this is the
   first such tab then tabwin->mbar will be NULL, otherwise
   if will be some other page's mbar, which will have to
   be swapped out
*/

void tabwin_register_toolbar (windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    gulong handler_id;

    /* take out a reference to @vwin's toolbar to prevent
       its auto-destruction; also ensure that the pointer
       goes to NULL on destruction
    */
    g_object_ref(vwin->mbar);
    handler_id = g_signal_connect(G_OBJECT(vwin->mbar), "destroy",
				  G_CALLBACK(gtk_widget_destroyed), 
				  &vwin->mbar);
    g_object_set_data(G_OBJECT(vwin->mbar), "destroy-id",
		      GUINT_TO_POINTER(handler_id));

#if TDEBUG
    fprintf(stderr, "*** register_toolbar: vwin=%p has toolbar=%p\n",
	    (void *) vwin, (void *) vwin->mbar);
#endif

    if (tabwin->mbar != NULL) {
	tabwin_remove_toolbar(tabwin);
    }

    tabwin_insert_toolbar(tabwin, vwin);
}

/* This is used, inter alia, for an "untitled" tab: 
   to set its real filename when it is saved 
*/

void tabwin_tab_set_title (windata_t *vwin, const char *title)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkWidget *tab, *label;

    tab = gtk_notebook_get_tab_label(GTK_NOTEBOOK(tabwin->tabs), 
				     vwin->main);

    label = g_object_get_data(G_OBJECT(tab), "label");
    if (label != NULL) {
	gtk_label_set_text(GTK_LABEL(label), title);
    }

    label = g_object_get_data(G_OBJECT(tab), "mlabel");
    if (label != NULL) {
	gtk_label_set_text(GTK_LABEL(label), title);
    }
}

/* set or unset the "modified flag" (trailing asterisk on
   the filename) for the tab label for a page in tabbed 
   editor
*/

void tabwin_tab_set_status (windata_t *vwin)
{
    gboolean unsaved = (vwin->flags & VWIN_CONTENT_CHANGED);
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkWidget *tab, *label;
    const gchar *text, *p;
    gchar *modtext = NULL;

    tab = gtk_notebook_get_tab_label(GTK_NOTEBOOK(tabwin->tabs), 
				     vwin->main);
    label = g_object_get_data(G_OBJECT(tab), "label");

    if (label != NULL) {
	text = gtk_label_get_text(GTK_LABEL(label));
	p = strstr(text, " *");
	if (unsaved && p == NULL) {
	    modtext = g_strdup_printf("%s *", text);
	} else if (!unsaved && p != NULL) {
	    modtext = g_strndup(text, strlen(text) - 2);
	}
    }

    if (modtext != NULL) {
	gtk_label_set_text(GTK_LABEL(label), modtext);
	label = g_object_get_data(G_OBJECT(tab), "mlabel");
	if (label != NULL) {
	    gtk_label_set_text(GTK_LABEL(label), modtext);
	}
	g_free(modtext);
    }
}

void show_tabbed_viewer (windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    int np = gtk_notebook_get_n_pages(notebook);

    gtk_widget_show_all(vwin->main);

    if (np > 1) {
	int pgnum = gtk_notebook_page_num(notebook, vwin->main);

	gtk_notebook_set_current_page(notebook, pgnum);
    }

#if GTK_MAJOR_VERSION == 2 && GTK_MAJOR_VERSION < 18
    if (!GTK_WIDGET_VISIBLE(tabwin->main)) {
	gtk_widget_show_all(tabwin->main);
    }
#else
    if (!gtk_widget_get_visible(tabwin->main)) {
	gtk_widget_show_all(tabwin->main);
    }
#endif

    if (vwin->role == EDIT_SCRIPT) {
	gtk_drag_dest_set(vwin->text,
			  GTK_DEST_DEFAULT_ALL,
			  tabwin_drag_targets, 1,
			  GDK_ACTION_COPY);
	g_signal_connect(G_OBJECT(vwin->text), "drag-data-received",
			 G_CALLBACK(tabwin_handle_drag),
			 tabwin);
    }

    gtk_window_present(GTK_WINDOW(tabwin->main));
}

/* move among the editor tabs via keyboard */

void tabwin_navigate (windata_t *vwin, guint key)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);

    if (key == GDK_less) {
	gtk_notebook_prev_page(notebook);
    } else if (key == GDK_greater) {
	gtk_notebook_next_page(notebook);
    } else {
	/* numeric value, 1 to 9 */
	gtk_notebook_set_current_page(notebook, key - 1);
    }
}

static void size_new_toplevel (windata_t *vwin)
{
    int cw = get_char_width(vwin->text);
    int hsize, vsize;

    if (vwin->role == EDIT_SCRIPT) {
	hsize = SCRIPT_WIDTH;
	vsize = SCRIPT_HEIGHT;
    } else {
	hsize = 63; /* MODEL_WIDTH ? */
	vsize = MODEL_HEIGHT;
    }

    hsize *= cw;
    hsize += 48;

    gtk_window_set_default_size(GTK_WINDOW(vwin->main), hsize, vsize);    
}

static gchar *title_from_vwin (windata_t *vwin)
{
    if (vwin->role == VIEW_MODEL) {
	MODEL *pmod = vwin->data;

	return g_strdup_printf(_("gretl: model %d"), pmod->ID);
    } else {
	return title_from_filename(vwin->fname, TRUE);
    }
}

/* show or hide the New and Open toolbar items, which occupy
   the first two slots on the toolbar
*/

void script_editor_show_new_open (windata_t *vwin, gboolean show)
{
    GtkToolItem *item0, *item1;

    item0 = gtk_toolbar_get_nth_item(GTK_TOOLBAR(vwin->mbar), 0);
    item1 = gtk_toolbar_get_nth_item(GTK_TOOLBAR(vwin->mbar), 1);

    if (show) {
	gtk_widget_show(GTK_WIDGET(item0));
	gtk_widget_show(GTK_WIDGET(item1));
    } else {
	gtk_widget_hide(GTK_WIDGET(item0));
	gtk_widget_hide(GTK_WIDGET(item1));
    }	
}

/* response to pulling a script or model out of the tabbed
   context: we need to give the content its own window
*/

void undock_tabbed_viewer (GtkWidget *w, windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    gint pg = gtk_notebook_page_num(notebook, vwin->main);
    GtkWidget *hbox = vwin->main;
    gulong handler_id;
    gchar *title;

    if (gtk_notebook_get_n_pages(notebook) < 2) {
	/* we won't do this if there's only one page in the
	   viewer */
	return;
    }

#if TDEBUG
    fprintf(stderr, "undock_tabbed_viewer: starting on vwin at %p\n",
	    (void *) vwin);
#endif

    /* disconnect stuff */
    vwin->main = vwin->topmain = NULL;

    /* remove signals and data from hbox (ex vwin->main) */
    g_object_steal_data(G_OBJECT(hbox), "vwin");
    g_signal_handlers_disconnect_by_func(hbox,
					 free_windata,
					 vwin);
    handler_id = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(vwin->mbar),
						    "destroy-id"));
    if (handler_id > 0) {
	g_signal_handler_disconnect(vwin->mbar, handler_id);
	g_object_steal_data(G_OBJECT(vwin->mbar), "destroy-id");
    }

    /* extract vwin->vbox from its tabbed holder */
    g_object_ref(vwin->vbox);
    gtk_container_remove(GTK_CONTAINER(hbox), vwin->vbox);
    gtk_notebook_remove_page(notebook, pg);

    /* tweak vbox params */
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);

    /* build new shell for @vwin */
    vwin->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    title = title_from_vwin(vwin);
    gtk_window_set_title(GTK_WINDOW(vwin->main), title);
    g_free(title);
    handler_id = g_signal_connect(G_OBJECT(vwin->main), "destroy", 
				  G_CALLBACK(free_windata), vwin);
    g_object_set_data(G_OBJECT(vwin->main), "destroy-id",
		      GUINT_TO_POINTER(handler_id));
    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
    size_new_toplevel(vwin);

#if TDEBUG
    fprintf(stderr, "*** undock_tabbed_viewer: new main=%p, mbar=%p\n", 
	    (void *) vwin->main, (void *) vwin->mbar);
#endif

    /* add box for toolbar, pack it, drop extra ref., then
       remove the "tabbed" flag (note that the tabbed flag
       is wanted so that vwin_pack_toolbar() will put the
       toolbar up top)
    */
    vwin_pack_toolbar(vwin);
    g_object_unref(vwin->mbar);
    vwin->flags &= ~VWIN_TABBED;

    if (vwin->role == EDIT_SCRIPT) {
	script_editor_show_new_open(vwin, FALSE);
    }

    /* put vbox into new top-level window and drop extra ref. */
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);
    g_object_unref(vwin->vbox);

    /* connect delete signal for single-script window */
    g_signal_connect(G_OBJECT(vwin->main), "delete-event", 
		     G_CALLBACK(query_save_text), vwin);

    /* add to window list and attach window-key signals */
    window_list_add(vwin->main, vwin->role);

    /* add key-catcher for single-item window */
    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_viewer_key), vwin);

#ifndef G_OS_WIN32
    set_wm_icon(vwin->main);
#endif

    gtk_widget_show(vwin->main);
    gtk_widget_grab_focus(vwin->text);

#if TDEBUG
    fprintf(stderr, "undock_tabbed_viewer: done\n");
#endif
}

static void dock_viewer (GtkWidget *w, windata_t *vwin)
{
    tabwin_t *tabwin;
    GtkWidget *oldmain;
    GtkWidget *box;
    gulong handler_id;
    gchar *info = NULL;

    tabwin = (vwin->role == VIEW_MODEL)? tabmod : tabedit;
    if (tabwin == NULL) {
	return;
    }

#if TDEBUG
    fprintf(stderr, "dock_viewer: starting on vwin at %p\n",
	    (void *) vwin);
#endif

    oldmain = vwin->main;
    gtk_widget_hide(oldmain);

    /* disconnect */
    vwin->main = NULL;

    /* remove data, and also remove destruction-related signals,
       from stand-alone vwin->main */
    g_object_steal_data(G_OBJECT(oldmain), "vwin");
    g_signal_handlers_disconnect_by_func(oldmain,
					 free_windata,
					 vwin);
    g_signal_handlers_disconnect_by_func(oldmain,
					 query_save_text,
					 vwin);

    /* grab info for title */
    if (vwin->role == EDIT_SCRIPT) {
	info = g_strdup(vwin->fname);
    } else {
	const gchar *tmp = gtk_window_get_title(GTK_WINDOW(oldmain));

	if (!strncmp(tmp, "gretl: ", 7)) {
	    tmp += 7;
	}
	info = g_strdup(tmp);
    }    

    /* extract vwin->vbox from oldmain and trash oldmain */
    g_object_ref(vwin->vbox);
    gtk_container_remove(GTK_CONTAINER(oldmain), vwin->vbox);
    gtk_widget_destroy(oldmain);

#if TDEBUG
    fprintf(stderr, "dock_viewer: vwin->vbox at %p\n", 
	    (void *) vwin->vbox);
#endif

    if (vwin->role == EDIT_SCRIPT) {
	script_editor_show_new_open(vwin, TRUE);
    }

    /* extract vwin->mbar */
    g_object_ref(vwin->mbar);
    box = gtk_widget_get_parent(vwin->mbar);
    gtk_container_remove(GTK_CONTAINER(box), vwin->mbar);

    /* create new vwin->main, etc. */
    vwin->main = gtk_hbox_new(FALSE, 0);
    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
    handler_id = g_signal_connect(G_OBJECT(vwin->main), "destroy", 
				  G_CALLBACK(free_windata), vwin);
    g_object_set_data(G_OBJECT(vwin->main), "destroy-id",
		      GUINT_TO_POINTER(handler_id));

    make_viewer_tab(tabwin, vwin, info);
    vwin->topmain = tabwin->main;
    vwin->flags = VWIN_TABBED;
    g_free(info);

    /* tweak vbox params and insert */
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 1);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 1);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);
    g_object_unref(vwin->vbox);

    /* repack toolbar in tabwin */
    tabwin_register_toolbar(vwin);
    g_object_unref(vwin->mbar);

    show_tabbed_viewer(vwin);
    gtk_widget_grab_focus(vwin->text);

#if TDEBUG
    fprintf(stderr, "dock_viewer: done\n");
#endif
}

gboolean window_is_undockable (windata_t *vwin)
{
    if (vwin->topmain != NULL) {
	tabwin_t *tabwin = vwin_get_tabwin(vwin);

	if (gtk_notebook_get_n_pages(GTK_NOTEBOOK(tabwin->tabs)) > 1) {
	    return TRUE;
	}
    }

    return FALSE;
}

gboolean window_is_dockable (windata_t *vwin)
{
    if (vwin->topmain == NULL) {
	if (vwin->role == EDIT_SCRIPT && tabedit != NULL) {
	    return TRUE;
	} else if (vwin->role == VIEW_MODEL && tabmod != NULL) {
	    return TRUE;
	}
    }

    return FALSE;
}

void add_undock_popup_item (GtkWidget *menu, windata_t *vwin)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(_("Move to new window"));
    g_signal_connect(G_OBJECT(item), "activate",
		     G_CALLBACK(undock_tabbed_viewer),
		     vwin);
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
}

void add_dock_popup_item (GtkWidget *menu, windata_t *vwin)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(_("Move to tabbed window"));
    g_signal_connect(G_OBJECT(item), "activate",
		     G_CALLBACK(dock_viewer),
		     vwin);
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
}

windata_t *tabwin_get_editor_for_file (const char *filename,
				       GtkWidget *w)
{
    windata_t *ret = NULL;

    if (w == tabedit->main) {
	GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
	int i, n = gtk_notebook_get_n_pages(notebook);
	GtkWidget *tab;
	windata_t *vwin;

	for (i=0; i<n; i++) {
	    tab = gtk_notebook_get_nth_page(notebook, i);
	    vwin = g_object_get_data(G_OBJECT(tab), "vwin");
	    if (vwin != NULL && !strcmp(filename, vwin->fname)) {
		ret = vwin;
		break;
	    }
	}
    }

    return ret;
}

void tabwin_tab_present (windata_t *vwin)
{
    tabwin_t *tabwin = g_object_get_data(G_OBJECT(vwin->topmain), 
					 "tabwin");
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    int i, n = gtk_notebook_get_n_pages(notebook);
    GtkWidget *tab;

    for (i=0; i<n; i++) {
	tab = gtk_notebook_get_nth_page(notebook, i);
	if (vwin == g_object_get_data(G_OBJECT(tab), "vwin")) {
	    gint pg = gtk_notebook_page_num(notebook, vwin->main);

	    gtk_notebook_set_current_page(notebook, pg);
	    break;
	}
    }

    gtk_window_present(GTK_WINDOW(vwin->topmain));
}

void tabwin_close_models_viewer (GtkWidget *w)
{
    if (tabmod != NULL && w == tabmod->main) {
	gtk_widget_destroy(w);
    }
}

static void tabwin_unregister_dialog (GtkWidget *w, tabwin_t *tabwin)
{
    if (tabwin != NULL && (tabwin == tabmod || tabwin == tabedit)) {
	/* @tabwin will be an invalid pointer if it
	   got a delete-event before execution gets
	   here
	*/
	GtkWidget *tab = tabwin->dlg_owner;

#if TDEBUG
	fprintf(stderr, "*** unregister_dialog: owner = %p\n", (void *) tab);
#endif

	if (tab != NULL) {
	    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
	    windata_t *vwin = g_object_get_data(G_OBJECT(tab), "vwin");

	    gtk_widget_set_sensitive(GTK_WIDGET(vwin->mbar), TRUE);
	    if (gtk_notebook_get_n_pages(notebook) > 1) {
		gtk_notebook_set_tab_detachable(notebook, tab, TRUE);
	    }
	    tabwin->dlg_owner = NULL;
	}
	tabwin->dialog = NULL;
    }
}

/* Called when a tabbed viewer spawns a dialog that becomes 
   invalid if the currently active tab is destroyed. We make 
   make the current tab undestroyable and undetachable for the 
   duration.
*/

void tabwin_register_dialog (GtkWidget *w, gpointer p)
{
    tabwin_t *tabwin = g_object_get_data(G_OBJECT(p), "tabwin");
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    gint pg = gtk_notebook_get_current_page(notebook);
    GtkWidget *tab = gtk_notebook_get_nth_page(notebook, pg);
    windata_t *vwin = g_object_get_data(G_OBJECT(tab), "vwin");

#if TDEBUG
    fprintf(stderr, "*** tabwin_register_dialog: w = %p, tab=%p\n", 
	    (void *) w, (void *) tab);
#endif

    gtk_widget_set_sensitive(vwin->mbar, FALSE);
    gtk_notebook_set_tab_detachable(notebook, tab, FALSE);

    tabwin->dialog = w;
    tabwin->dlg_owner = tab;

    g_signal_connect(G_OBJECT(w), "destroy",
		     G_CALLBACK(tabwin_unregister_dialog), 
		     tabwin);
}

int viewer_n_siblings (windata_t *vwin) 
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    int n = 0;

    if (tabwin != NULL) {
	n = gtk_notebook_get_n_pages(GTK_NOTEBOOK(tabwin->tabs));
	if (n > 0) n--;
    }

    return n;
}

int highest_numbered_var_in_tabwin (tabwin_t *tabwin, 
				    const DATASET *dset)
{
    int vmax = 0;

    if (tabwin->role == VIEW_MODEL) {
	GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
	int n = gtk_notebook_get_n_pages(notebook);
	const MODEL *pmod;
	windata_t *vwin;
	GtkWidget *tab;
	int i, m_vmax;
	
	for (i=0; i<n; i++) {
	    tab = gtk_notebook_get_nth_page(notebook, i);
	    vwin = g_object_get_data(G_OBJECT(tab), "vwin");
	    pmod = vwin->data;
	    m_vmax = highest_numbered_var_in_model(pmod, dset);
	    if (m_vmax > vmax) {
		vmax = m_vmax;
	    }
	}
    }

    return vmax;
}

windata_t *window_get_active_vwin (GtkWidget *window) 
{
    windata_t *vwin = g_object_get_data(G_OBJECT(window), "vwin");

    if (vwin == NULL) {
	tabwin_t *tabwin = g_object_get_data(G_OBJECT(window), "tabwin");

	if (tabwin != NULL) {
	    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
	    gint pg = gtk_notebook_get_current_page(notebook);
	    GtkWidget *tab = gtk_notebook_get_nth_page(notebook, pg);

	    vwin = g_object_get_data(G_OBJECT(tab), "vwin");
	}
    }

    return vwin;
}
