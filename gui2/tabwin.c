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

#define TDEBUG 0

typedef struct tabwin_t_  tabwin_t;

struct tabwin_t_ {
    int role;         /* what's tabwin doing? */ 
    GtkWidget *main;  /* top-level GTK window */
    GtkWidget *mbox;  /* horizontal box to hold menu bar */
    GtkWidget *mbar;  /* menu bar */
    GtkWidget *tabs;  /* notebook for tabs */
};

/* We support one tabbed editor, for gretl scripts --
   and we may support a tabbed viewer for models */

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
    GtkWidget *button = g_object_get_data(G_OBJECT(lbl), "button");

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

    if (np < 5) {
	gtk_notebook_popup_disable(notebook);
    }    

    if (np == 1) {
	/* only one tab left after removal: this page should
	   not display its own closer button, nor should it
	   be detachable.
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

static void editor_tab_destroy (GtkWidget *w, windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    gint pg = gtk_notebook_page_num(notebook, vwin->main);

#if TDEBUG
    fprintf(stderr, "*** editor_tab_destroy: vwin = %p\n", (void *) vwin);
#endif

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
    g_object_unref(G_OBJECT(vwin->mbar));

    gtk_notebook_remove_page(notebook, pg);
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
    g_signal_connect(button, "clicked", G_CALLBACK(editor_tab_destroy), 
		     vwin);
    gtk_container_add(GTK_CONTAINER(tab), button);
    g_object_set_data(G_OBJECT(tab), "button", button);
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
	return g_strdup_printf("untitled(%d)", idx);
    } else {
	return g_strdup("untitled");
    }
}

/* create and add tab with filename and closer button */

static GtkWidget *make_viewer_tab (tabwin_t *tabwin, 
				   windata_t *vwin, 
				   const gchar *info,
				   int starting)
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
    gtk_widget_show_all(tab);

    notebook = GTK_NOTEBOOK(tabwin->tabs);
    gtk_notebook_append_page_menu(notebook, vwin->main, 
				  tab, mlabel);

    return tab;
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
    g_signal_connect(tabwin->tabs, "switch-page",
		     G_CALLBACK(switch_page_callback), tabwin);
    g_signal_connect(tabwin->tabs, "create-window",
		     G_CALLBACK(detach_tab_callback), tabwin);
    g_signal_connect(tabwin->tabs, "page-added",
		     G_CALLBACK(page_added_callback), tabwin);
    g_signal_connect(tabwin->tabs, "page-removed",
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
	if (tabedit == NULL) {
	    *starting = 1;
	    tabedit = tabwin = make_tabbed_viewer(role);
	} else {
	    tabwin = tabedit;
	}
    } else if (role == VIEW_MODEL) {
	if (tabmod == NULL) {
	    *starting = 1;
	    tabmod = tabwin = make_tabbed_viewer(role);
	} else {
	    tabwin = tabmod;
	}
    }	

    return tabwin;
}

/* Create an editor tab, as an alternative to a stand-alone editor
   window. We build the tabbed editor if need be, otherwise we
   stick a new tab into the existing editor.
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
    g_object_set_data(G_OBJECT(vwin->main), "handler-id",
		      GUINT_TO_POINTER(handler_id));

#if TDEBUG
    fprintf(stderr, "*** viewer_tab_new: vwin=%p, main hbox=%p\n", 
	    (void *) vwin, (void *) vwin->main);
#endif

    make_viewer_tab(tabwin, vwin, info, starting);
    vwin->topmain = tabwin->main;

    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 1);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    if (starting) {
	add_window_list_item(tabwin->main, role);
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

    /* take out a reference to @vwin's toolbar to prevent
       its auto-destruction; also ensure that the pointer
       goes to NULL on destruction
    */
    g_object_ref(G_OBJECT(vwin->mbar));
    g_signal_connect(G_OBJECT(vwin->mbar), "destroy",
                     G_CALLBACK(gtk_widget_destroyed), &vwin->mbar);

#if TDEBUG
    fprintf(stderr, "*** register_toolbar: vwin=%p has toolbar=%p\n",
	    (void *) vwin, (void *) vwin->mbar);
#endif

    if (tabwin->mbar != NULL) {
	tabwin_remove_toolbar(tabwin);
    }

    tabwin_insert_toolbar(tabwin, vwin);
}

/* for an "untitled" tab: set its real filename when it is
   saved */

void tabwin_set_tab_title (windata_t *vwin, gchar *fname)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkWidget *tab, *label;

    tab = gtk_notebook_get_tab_label(GTK_NOTEBOOK(tabwin->tabs), 
				     vwin->main);

    label = g_object_get_data(G_OBJECT(tab), "label");
    if (label != NULL) {
	gtk_label_set_text(GTK_LABEL(label), fname);
    }

    label = g_object_get_data(G_OBJECT(tab), "mlabel");
    if (label != NULL) {
	gtk_label_set_text(GTK_LABEL(label), fname);
    }
}

/* set or unset the "modified flag" (trailing asterisk on
   the filename) for the tab label for a page in tabbed 
   editor
*/

void tabwin_set_tab_status (windata_t *vwin)
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

static void size_new_toplevel (GtkWidget *w, windata_t *vwin)
{
    int cw = get_char_width(vwin->text);
    int hsize = 78 * cw + 48;
    int vsize = 370;

    gtk_window_set_default_size(GTK_WINDOW(w), hsize, vsize);    
}

void undock_tabbed_viewer (GtkWidget *w, windata_t *vwin)
{
    tabwin_t *tabwin = vwin_get_tabwin(vwin);
    GtkNotebook *notebook = GTK_NOTEBOOK(tabwin->tabs);
    gint pg = gtk_notebook_page_num(notebook, vwin->main);
    gulong handler_id;
    GtkWidget *mainwin;
    gchar *title;

    /* we'll not do this if there's only one page in the
       editor
    */
    if (gtk_notebook_get_n_pages(notebook) < 2) {
	return;
    }

    /* take a reference to the guts of @vwin, undo the attachment
       to its holder, and extract from the tabbed context
    */
    g_object_ref(G_OBJECT(vwin->vbox));
    g_object_set_data(G_OBJECT(vwin->main), "vwin", NULL);
    handler_id = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(vwin->main),
						    "handler-id"));
    g_signal_handler_disconnect(vwin->main, handler_id);
    gtk_container_remove(GTK_CONTAINER(vwin->main), vwin->vbox);
    gtk_notebook_remove_page(notebook, pg);

    /* build new shell for @vwin's vbox */
    mainwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    title = title_from_filename(vwin->fname, TRUE);
    gtk_window_set_title(GTK_WINDOW(mainwin), title);
    g_free(title);
    g_signal_connect(G_OBJECT(mainwin), "destroy", 
		     G_CALLBACK(free_windata), vwin);
    g_object_set_data(G_OBJECT(mainwin), "vwin", vwin);
    size_new_toplevel(mainwin, vwin);

    /* disconnect tabbed top-level */
    vwin->topmain = NULL;

    /* tweak vbox params */
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);

    /* add box for toolbar, pack it, drop extra ref., then
       remove the "tabbed" flag (note that the tabbed flag
       is wanted so that vwin_pack_toolbar() will put the
       toolbar up top)
    */
    vwin_pack_toolbar(vwin);
    g_object_unref(G_OBJECT(vwin->mbar));
    vwin->flags &= ~VWIN_TABBED;

    /* put vbox into new top-level window and drop extra ref. */
    gtk_container_add(GTK_CONTAINER(mainwin), vwin->vbox);
    g_object_unref(G_OBJECT(vwin->vbox));

     /* connect toplevel and set signals */
    vwin->main = mainwin;
    g_signal_connect(G_OBJECT(vwin->main), "delete-event", 
		     G_CALLBACK(query_save_text), vwin);
    g_object_set_data(G_OBJECT(vwin->main), "role", 
		      GINT_TO_POINTER(vwin->role));
    winstack_add(vwin->main);

    add_window_list_item(vwin->main, vwin->role);

#ifndef G_OS_WIN32
    set_wm_icon(vwin->main);
#endif

    gtk_widget_show_all(vwin->main);
    gtk_widget_grab_focus(vwin->text);
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
