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
    GtkWidget *mbar;      /* menu bar */
    GtkWidget *tabs;      /* notebook for tabs */
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
    if (tabedit->mbar == vwin->mbar) {
	tabedit->mbar = NULL;
    }
    gtk_widget_destroy(vwin->mbar);
    vwin->mbar = NULL;

    gtk_notebook_remove_page(notebook, pg);

    if (np == 2) {
	/* so only one page left after removal */
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
	title = g_strdup("untitled");
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
   toolbar into place in tabedit (and hide the old
   one)
*/

static gboolean tabedit_switch_page (GtkNotebook *tabs,
				     GtkNotebook *page,
				     gint pgnum,
				     tabwin_t *twin)
{
    if (page != NULL) {
	windata_t *vwin;

	vwin = g_object_get_data(G_OBJECT(page), "vwin");
	if (vwin != NULL) {
	    if (tabedit->mbar != NULL) {
		gtk_widget_hide(tabedit->mbar);
	    }
	    tabedit->mbar = vwin->mbar;
	    if (tabedit->mbar != NULL) {
		gtk_widget_show(tabedit->mbar);
	    }
	}
    }

    return FALSE;
}

/* build the tabbed editor for gretl scripts */

static tabwin_t *make_tabedit (void)
{
    GtkWidget *hbox, *vbox;

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

    /* hbox to hold menu bar */
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    g_object_set_data(G_OBJECT(tabedit->main), "hbox", hbox);
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
    g_signal_connect(G_OBJECT(vwin->main), "destroy", 
		     G_CALLBACK(free_windata), vwin);

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

void tabwin_register_toolbar (windata_t *vwin)
{
    GtkWidget *hbox;

    hbox = g_object_get_data(G_OBJECT(vwin->topmain), "hbox");
    if (tabedit->mbar != NULL) {
	gtk_widget_hide(tabedit->mbar);
    }
    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, TRUE, TRUE, 0);
    tabedit->mbar = vwin->mbar;
    gtk_widget_show_all(tabedit->mbar);
}

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
