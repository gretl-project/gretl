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

typedef struct tabwin_t_  tabwin_t;

struct tabwin_t_ {
    GtkWidget *main;      /* top-level GTK window */
    GtkWidget *vbox;      /* vertical container */
    GtkWidget *mbar;      /* menu bar */
    GtkWidget *tabs;      /* notebook for tabs */
};

/* We only support one tabbed editor, for gretl scripts */

static tabwin_t *tabedit;

static void tabedit_destroy (GtkWidget *w, gpointer data)
{
    /* FIXME more stuff needed here? */
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

void maybe_destroy_tabwin (GtkWidget *tmain)
{
    if (!maybe_block_tabedit_quit()) {
	gtk_widget_destroy(tmain);
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

/* active the tab's closer button */

static void viewer_tab_show_closer (GtkWidget *tab)
{
    GtkWidget *button = g_object_get_data(G_OBJECT(tab), "button");

    if (button != NULL) {
	gtk_widget_show(button);
    }
}

/* callback for tab-specific close button */

static void editor_tab_destroy (GtkWidget *w, GtkWidget *vmain)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
    gint pg = gtk_notebook_page_num(notebook, vmain);
    gint np = gtk_notebook_get_n_pages(notebook);

    if (np == 1) {
	/* shouldn't happen */
	gtk_widget_destroy(tabedit->main);
    } else {
	gtk_notebook_remove_page(notebook, pg);
	if (np == 2) {
	    GtkWidget *child = gtk_notebook_get_nth_page(notebook, 0);
	    GtkWidget *tab = gtk_notebook_get_tab_label(notebook, child);

	    viewer_tab_hide_closer(tab);
	}
    }
}

/* put a tab-specific close button next to the tab's label */

static void viewer_tab_add_closer (GtkWidget *tab, GtkWidget *vmain)
{
    GtkWidget *img, *button;
	
    img = gtk_image_new_from_stock(GTK_STOCK_CLOSE, 
				   GTK_ICON_SIZE_MENU);
    button = gtk_button_new();
    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    gtk_container_set_border_width(GTK_CONTAINER(button), 0);
    gtk_container_add(GTK_CONTAINER(button), img);
    g_signal_connect(button, "clicked", G_CALLBACK(editor_tab_destroy), 
		     vmain);
    gtk_container_add(GTK_CONTAINER(tab), button);
    g_object_set_data(G_OBJECT(tab), "button", button);
}

static GtkWidget *make_viewer_tab (tabwin_t *twin, 
				   windata_t *vwin, 
				   const gchar *filename,
				   int starting)
{
    gchar *title = NULL;
    const gchar *p;
    GtkWidget *tab;
    GtkWidget *label;

    tab = gtk_hbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(tab), 0);

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

    viewer_tab_add_closer(tab, vwin->main);
    gtk_widget_show_all(tab);

    if (starting) {
	viewer_tab_hide_closer(tab);
    }
    
    return tab;
}

static void tab_box_config (windata_t *vwin)
{
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 1);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);
    
#ifndef G_OS_WIN32
    set_wm_icon(vwin->topmain);
#endif
}

#if 0
static void retarget_toolbar (void)
{
    GtkToolbar *tbar = GTK_TOOLBAR(tabedit->tbar);
    GtkToolItem *item;
    int i, n;
    
    n = gtk_toolbar_get_n_items(tbar);

    for (i=0; i<n; i++) {
	item = gtk_toolbar_get_nth_item(tbar, i);
	g_signal_connect(item, "clicked", func, data);
    }
}
#endif

/* on switching the current page, connect the new vwin
   to tabedit's toolbar 
*/

static gboolean tabedit_switch_page (GtkNotebook *tabs,
				     gpointer arg1,
				     gint newpage,
				     tabwin_t *twin)
{
    GtkWidget *w = gtk_notebook_get_nth_page(tabs, newpage);
    windata_t *vwin;

    if (w != NULL) {
	vwin = g_object_get_data(G_OBJECT(w), "vwin");
	if (vwin != NULL) {
	    vwin->mbar = tabedit->mbar;
	}
    }

    return FALSE;
}

/* build the tabbed editor for gretl scripts */

static tabwin_t *make_tabedit (void)
{
    GtkWidget *hbox;

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
    tabedit->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(tabedit->vbox), 0);
    gtk_container_set_border_width(GTK_CONTAINER(tabedit->vbox), 0);
    gtk_container_add(GTK_CONTAINER(tabedit->main), tabedit->vbox);

    /* hbox to hold menu bar */
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tabedit->vbox), hbox, FALSE, FALSE, 0);
    g_object_set_data(G_OBJECT(tabedit->main), "hbox", hbox);
    tabedit->mbar = NULL;

    /* notebook */
    tabedit->tabs = gtk_notebook_new();
    g_signal_connect(tabedit->tabs, "switch-page",
		     G_CALLBACK(tabedit_switch_page), tabedit);
    gtk_container_add(GTK_CONTAINER(tabedit->vbox), tabedit->tabs);

    return tabedit;
}

/* Create an editor tab, as an alternative to a stand-alone editor
   window. We build the tabed editor if need be, otherwise we
   stick a new tab into the existing editor.
*/

windata_t *editor_tab_new (const char *filename)
{
    GtkNotebook *notebook;
    windata_t *vwin;
    GtkWidget *tab;
    int starting = 0;

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

    vwin->main = gtk_hbox_new(FALSE, 0);
    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
    g_signal_connect(G_OBJECT(vwin->main), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    tab = make_viewer_tab(tabedit, vwin, filename, starting);
    notebook = GTK_NOTEBOOK(tabedit->tabs);
    gtk_notebook_append_page(notebook, vwin->main, tab);
    vwin->topmain = tabedit->main;

    tab_box_config(vwin);

    if (starting) {
	add_window_list_item(tabedit->main, EDIT_SCRIPT); /* FIXME */
    } 

    return vwin;
}

void tabwin_register_toolbar (windata_t *vwin)
{
    GtkWidget *hbox;

    hbox = g_object_get_data(G_OBJECT(vwin->topmain), "hbox");
    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, TRUE, TRUE, 0);
    tabedit->mbar = vwin->mbar;
}

void show_tabbed_viewer (GtkWidget *vmain)
{
    GtkNotebook *notebook = GTK_NOTEBOOK(tabedit->tabs);
    int np = gtk_notebook_get_n_pages(notebook);

    gtk_widget_show_all(vmain);

    if (np > 1) {
	int pg = gtk_notebook_page_num(notebook, vmain);

	gtk_notebook_set_current_page(notebook, pg);

	if (pg == 1) {
	    GtkWidget *child = gtk_notebook_get_nth_page(notebook, 0);
	    GtkWidget *tab = gtk_notebook_get_tab_label(notebook, child);

	    viewer_tab_show_closer(tab);
	}
    }

    if (!gtk_widget_get_visible(tabedit->main)) {
	gtk_widget_show_all(tabedit->main);
    }
}
