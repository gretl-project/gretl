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

#define TDEBUG 0

typedef struct tabwin_t_  tabwin_t;

struct tabwin_t_ {
    GtkWidget *main;  /* top-level GTK window */
    GtkWidget *mbox;  /* hbox to hold menu bar */
    GtkWidget *mbar;  /* menu bar */
    GtkWidget *tabs;  /* notebook for tabs */
};

/* We only support one tabbed editor, for gretl scripts */

static tabwin_t *tabedit;

static void tabedit_destroy (GtkWidget *w, gpointer data)
{
    free(tabedit);
    tabedit = NULL;
}

static gboolean maybe_block_tabedit_quit (void)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
    int np = gtk_notebook_get_n_pages(notebook);
    gboolean ret = FALSE;

    if (np > 1) {
	gchar *msg = g_strdup_printf(_("Editing %d scripts: really quit?"), np);
	gint resp = yes_no_dialog(_("gretl: script editor"), msg, 0);

	if (resp != GRETL_YES) {
	    ret = TRUE;
	}
	g_free(msg);
    } else if (np == 1) {
	GtkWidget *page = gtk_notebook_get_nth_page(notebook, 0);

	if (page != NULL) {
	    windata_t *vwin = g_object_get_data(G_OBJECT(page), "vwin");
    
	    if (vwin_content_changed(vwin)) {
		ret = query_save_text(NULL, NULL, vwin);
	    }
	}
    }

    return ret;
}

static gboolean tabedit_quit_check (GtkWidget *w, GdkEvent *event, 
				    tabwin_t *twin)
{
    return maybe_block_tabedit_quit();
}

void maybe_destroy_tabwin (windata_t *vwin)
{
    if (!maybe_block_tabedit_quit()) {
	gtk_widget_destroy(vwin->topmain);
    }
}

/* called if there's only one tab left: don't show a tab-specific
   close button */

static void viewer_tab_hide_closer (GtkWidget *tab)
{
    GtkWidget *button = g_object_get_data(G_OBJECT(tab), "button");

    if (button != NULL) {
	gtk_widget_hide(button);
    }
}

/* activate the tab's closer button */

static void viewer_tab_show_closer (GtkWidget *tab)
{
    GtkWidget *button = g_object_get_data(G_OBJECT(tab), "button");

    if (button != NULL) {
	gtk_widget_show(button);
    }
}

static void tabedit_remove_toolbar (void)
{
#if TDEBUG
    fprintf(stderr, "*** removing toolbar at %p\n", (void *) tabedit->mbar);
#endif
    gtk_container_remove(GTK_CONTAINER(tabedit->mbox), 
			 tabedit->mbar);
    tabedit->mbar = NULL;
}

static void tabedit_insert_toolbar (windata_t *vwin)
{
#if TDEBUG
    fprintf(stderr, "*** inserting toolbar at %p\n", (void *) vwin->mbar);
#endif
    gtk_box_pack_start(GTK_BOX(tabedit->mbox), vwin->mbar, 
		       TRUE, TRUE, 0);
    gtk_widget_show_all(vwin->mbar);
    tabedit->mbar = vwin->mbar;
}

/* callback for tab-specific close button: this should be
   invoked only if there's at least one tab left after
   trashing the selected one
*/

static void editor_tab_destroy (GtkWidget *w, windata_t *vwin)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
    gint pg = gtk_notebook_page_num(notebook, vwin->main);
    gint np = gtk_notebook_get_n_pages(notebook);

    /* note: vwin->mbar is packed under tabedit, so it will not
       get destroyed automatically when the page is removed
    */
    if (tabedit->mbar != NULL && tabedit->mbar == vwin->mbar) {
	tabedit_remove_toolbar();
    }

    /* relinquish the extra reference */
#if TDEBUG
    fprintf(stderr, "*** unreffing toolbar at %p\n", (void *) vwin->mbar);
#endif
    g_object_unref(vwin->mbar);

    gtk_notebook_remove_page(notebook, pg);

    if (np == 2) {
	/* so only one page left after removal: this page should
	   not display its own closer button 
	*/
	GtkWidget *page = gtk_notebook_get_nth_page(notebook, 0);
	GtkWidget *tab = gtk_notebook_get_tab_label(notebook, page);

	viewer_tab_hide_closer(tab);
    }

    if (np < 5) {
	gtk_notebook_popup_disable(notebook);
    }
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

static gchar *untitled_title (tabwin_t *twin)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(twin->tabs);
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

/* add tab with filename and closer button */

static GtkWidget *make_viewer_tab (tabwin_t *twin, 
				   windata_t *vwin, 
				   const gchar *filename,
				   int starting)
{
    gchar *title = NULL;
    const gchar *p;
    GtkWidget *tab;
    GtkWidget *label;
    GtkWidget *mlabel;

    tab = gtk_hbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(tab), 0);

    if (strstr(filename, "script_tmp") != NULL) {
	title = untitled_title(twin);
    } else if ((p = strrchr(filename, SLASH)) != NULL) {
	title = g_strdup(p + 1);
    } else {
	title = g_strdup(filename);
    }

    label = gtk_label_new(title);
    mlabel = gtk_label_new(title);
    gtk_container_add(GTK_CONTAINER(tab), label);
    g_object_set_data(G_OBJECT(tab), "label", label);
    g_object_set_data(G_OBJECT(tab), "mlabel", mlabel);
    g_free(title);

    viewer_tab_add_closer(tab, vwin);
    gtk_widget_show_all(tab);

    if (starting) {
	viewer_tab_hide_closer(tab);
    }

    gtk_notebook_append_page_menu(GTK_NOTEBOOK(twin->tabs), 
				  vwin->main, tab, mlabel);

    return tab;
}

/* on switching the current page, put the new page's
   toolbar into place in tabedit (and remove the old
   one, if present)
*/

static gboolean tabedit_switch_page (GtkNotebook *tabs,
				     GtkNotebook *page,
				     gint pgnum,
				     tabwin_t *twin)
{
    windata_t *vwin = NULL;

    if (page != NULL) {
	vwin = g_object_get_data(G_OBJECT(page), "vwin");
    }

    if (vwin != NULL) {
	if (tabedit->mbar != NULL && tabedit->mbar != vwin->mbar) {
	    /* there's an "old" toolbar in place */
	    tabedit_remove_toolbar();
	}
	if (vwin->mbar != NULL && vwin->mbar != tabedit->mbar) {
	    /* a "new" toolbar should be shown */
	    tabedit_insert_toolbar(vwin);
	}
    }    

    return FALSE;
}

/* build the tabbed editor for gretl scripts */

static tabwin_t *make_tabedit (void)
{
    GtkWidget *vbox;

    tabedit = mymalloc(sizeof *tabedit);

    if (tabedit == NULL) {
	return NULL;
    }

    /* top-level window */
    tabedit->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(tabedit->main), 
			 _("gretl: script editor"));
    g_signal_connect(G_OBJECT(tabedit->main), "delete-event",
		     G_CALLBACK(tabedit_quit_check), tabedit);
    g_signal_connect(G_OBJECT(tabedit->main), "destroy", 
		     G_CALLBACK(tabedit_destroy), tabedit);
    g_object_set_data(G_OBJECT(tabedit->main), "tabedit", tabedit);

    /* vertically oriented container */
    vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vbox), 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 0);
    gtk_container_add(GTK_CONTAINER(tabedit->main), vbox);

    /* box to hold menu bar */
    tabedit->mbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), tabedit->mbox, FALSE, FALSE, 0);
    tabedit->mbar = NULL;

    /* notebook */
    tabedit->tabs = gtk_notebook_new();
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(tabedit->tabs), TRUE);
    g_signal_connect(tabedit->tabs, "switch-page",
		     G_CALLBACK(tabedit_switch_page), tabedit);
    gtk_container_add(GTK_CONTAINER(vbox), tabedit->tabs);

#ifndef G_OS_WIN32
    set_wm_icon(tabedit->main);
#endif

    return tabedit;
}

/* Create an editor tab, as an alternative to a stand-alone editor
   window. We build the tabed editor if need be, otherwise we
   stick a new tab into the existing editor.
*/

windata_t *editor_tab_new (const char *filename)
{
    windata_t *vwin;
    int starting = 0;
    gulong handler_id;

    if (tabedit == NULL) {
	starting = 1;
	tabedit = make_tabedit();
	if (tabedit == NULL) {
	    return NULL;
	}
    }

    vwin = vwin_new(EDIT_SCRIPT, NULL);
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

    make_viewer_tab(tabedit, vwin, filename, starting);
    vwin->topmain = tabedit->main;

    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 1);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    if (starting) {
	add_window_list_item(tabedit->main, EDIT_SCRIPT);
    } else {
	GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);

	if (gtk_notebook_get_n_pages(notebook) > 5) {
	    gtk_notebook_popup_enable(notebook);
	}
    }

    return vwin;
}

/* called when a new editor tab is added: if this is the
   first such tab then tabedit->mbar will be NULL, otherwise
   if will be some other page's mbar, which will have to
   be swapped out
*/

void tabwin_register_toolbar (windata_t *vwin)
{
    /* take out a reference to @vwin's toolbar to prevent
       its auto-destruction; also ensure that the pointer
       goes to NULL on destruction
    */
    g_object_ref(G_OBJECT(vwin->mbar));
    g_signal_connect(G_OBJECT(vwin->mbar), "destroy",
                     G_CALLBACK(gtk_widget_destroyed), &vwin->mbar);

#if TDEBUG
    fprintf(stderr, "*** vwin at %p: has toolbar at %p\n",
	    (void *) vwin, (void *) vwin->mbar);
#endif

    if (tabedit->mbar != NULL) {
	tabedit_remove_toolbar();
    }

    tabedit_insert_toolbar(vwin);
}

/* for an "untitled" tab: set its real filename when it is
   saved */

void tabwin_set_tab_title (windata_t *vwin, gchar *fname)
{
    GtkWidget *tab, *label;

    tab = gtk_notebook_get_tab_label(GTK_NOTEBOOK(tabedit->tabs), 
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
    GtkWidget *tab, *label;
    const gchar *text, *p;
    gchar *modtext = NULL;

    tab = gtk_notebook_get_tab_label(GTK_NOTEBOOK(tabedit->tabs), 
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

void show_tabbed_viewer (GtkWidget *vmain)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
    int np = gtk_notebook_get_n_pages(notebook);

    gtk_widget_show_all(vmain);

    if (np > 1) {
	int pgnum = gtk_notebook_page_num(notebook, vmain);

	gtk_notebook_set_current_page(notebook, pgnum);

	if (pgnum == 1) {
	    GtkWidget *page = gtk_notebook_get_nth_page(notebook, 0);
	    GtkWidget *tab = gtk_notebook_get_tab_label(notebook, page);

	    viewer_tab_show_closer(tab);
	}
    }

#if GTK_MAJOR_VERSION == 2 && GTK_MAJOR_VERSION < 18
    if (!GTK_WIDGET_VISIBLE(tabedit->main)) {
	gtk_widget_show_all(tabedit->main);
    }
#else
    if (!gtk_widget_get_visible(tabedit->main)) {
	gtk_widget_show_all(tabedit->main);
    }
#endif
}

/* move left or right among the editor tabs via keyboard */

void tabwin_navigate (windata_t *vwin, guint key)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
    int pgnum = gtk_notebook_get_current_page(notebook);

    if (key == GDK_less) {
	if (pgnum > 0) {
	    gtk_notebook_set_current_page(notebook, pgnum - 1);
	}
    } else if (key == GDK_greater) {
	int np = gtk_notebook_get_n_pages(notebook);

	if (pgnum < np - 1) {
	    gtk_notebook_set_current_page(notebook, pgnum + 1);
	}
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
    GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
    gint pg = gtk_notebook_page_num(notebook, vwin->main);
    gulong handler_id;
    GtkWidget *mainwin;

    /* we'll not do this is there's only one page in the
       editor
    */
    if (gtk_notebook_get_n_pages(notebook) < 2) {
	return;
    }

    /* take a reference to the guts of @vwin, undo the attachment
       to its holder, and extract from the tabbed context
    */
    g_object_ref(vwin->vbox);
    g_object_set_data(G_OBJECT(vwin->main), "vwin", NULL);
    handler_id = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(vwin->main),
						    "handler-id"));
    g_signal_handler_disconnect(vwin->main, handler_id);
    gtk_container_remove(GTK_CONTAINER(vwin->main), vwin->vbox);
    gtk_notebook_remove_page(notebook, pg);

    /* build new shell for @vwin's vbox */
    mainwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(mainwin), "gretl: oof!");
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
    g_object_unref(vwin->mbar);
    vwin->flags &= ~VWIN_TABBED;

    /* put vbox into new top-level window and drop extra ref. */
    gtk_container_add(GTK_CONTAINER(mainwin), vwin->vbox);
    g_object_unref(vwin->vbox);

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
	if (gtk_notebook_get_n_pages(GTK_NOTEBOOK(tabedit->tabs)) > 1) {
	    return TRUE;
	}
    }

    return FALSE;
}
