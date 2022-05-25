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

/* editbar.c: cut-dowm variant of toolbar.c for editor mode */

#include "gretl.h"
#include "textbuf.h"
#include "textutil.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "winstack.h"
#include "tabwin.h"
#include "gretl_www.h"
#include "editbar.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#include <gio/gio.h>

/* for viewer window toolbars */
#include "../pixmaps/mini.en.xpm"

/* for pop-up search entry */
#include "../pixmaps/close_16.xpm"

/* for window-finder menu */
#include "../pixmaps/mini.gretl.xpm"
#include "../pixmaps/mini.table.xpm"
#include "../pixmaps/mini.page.xpm"

enum {
    SAVE_ITEM = 1,
    SAVE_AS_ITEM,
    EDIT_ITEM,
    EXEC_ITEM,
    COPY_ITEM,
    PRINT_ITEM,
    HELP_ITEM,
    CMD_HELP_ITEM,
    EDIT_HANSL_ITEM,
    OPEN_ITEM,
    SPLIT_H_ITEM,
    SPLIT_V_ITEM,
    EDITOR_ITEM,
    NEW_ITEM,
    FIND_ITEM,
    CLOSE_ITEM
} viewbar_flags;

int toolbar_icon_size = GTK_ICON_SIZE_MENU;

struct png_stock_maker {
    char *fname;
    const char *id;
    gint8 in_menu;
};

struct png_stock_maker png_stocks[] = {
    { "winlist.png",   GRETL_STOCK_WINLIST },
    { "tools.png",     GRETL_STOCK_TOOLS, 0 },
    { "pushpin.png",   GRETL_STOCK_PIN, 0 },
    { "split_h.png",   GRETL_STOCK_SPLIT_H, 0 },
    { "join_h.png",    GRETL_STOCK_JOIN_H, 0 },
    { "split_v.png",   GRETL_STOCK_SPLIT_V, 0 },
    { "join_v.png",    GRETL_STOCK_JOIN_V, 0 },
    { "query.png",     GRETL_STOCK_QUERY, 0 },
};

struct xpm_stock_maker {
    char **xpm;
    const char *id;
};

#if GTK_MAJOR_VERSION == 3

static void try_auto_icon_sizing (int *bigger)
{
    GdkDisplay *display = gdk_display_get_default();
    GdkMonitor *monitor = gdk_display_get_primary_monitor(display);
    GdkScreen *screen = gdk_screen_get_default();

    if (monitor != NULL) {
	int mmw = gdk_monitor_get_width_mm(monitor);
	int pxw = gdk_screen_get_width(screen);
	double mm16 = 16 * mmw / (double) pxw;

	if (mm16 < 2.8) {
	    /* size of 16 pixels in millimeters */
	    fprintf(stderr, " auto-setting larger icons\n");
	    *bigger = 1;
	}
    }
}

#endif /* GTK3 */

void gretl_stock_icons_init (void)
{
    struct xpm_stock_maker xpm_stocks[] = {
	{ mini_en_xpm, GRETL_STOCK_EN },
	{ mini_gretl_xpm, GRETL_STOCK_GRETL},
	{ mini_table_xpm, GRETL_STOCK_TABLE},
	{ mini_page_xpm, GRETL_STOCK_PAGE},
	{ close_16_xpm, GRETL_STOCK_CLOSE}
    };
    static GtkIconFactory *gretl_factory;
    int n1 = G_N_ELEMENTS(png_stocks);
    int n2 = G_N_ELEMENTS(xpm_stocks);

    if (gretl_factory == NULL) {
	int bigger = (get_icon_sizing() == ICON_SIZE_MEDIUM);
	char icon_path[48], menu_path[48];
	gchar *p, *pm, *respath;
	GResource *icons;
	GtkIconSource *isrc;
	GtkIconSet *iset;
	GdkPixbuf *pbuf;
	int i;

#if GTK_MAJOR_VERSION == 3
	if (get_icon_sizing() == ICON_SIZE_AUTO) {
	    try_auto_icon_sizing(&bigger);
	}
#endif

	if (bigger) {
#if GTK_MAJOR_VERSION == 3
	    toolbar_icon_size = GTK_ICON_SIZE_LARGE_TOOLBAR;
#else
	    toolbar_icon_size = GTK_ICON_SIZE_SMALL_TOOLBAR;
#endif
	}

	gretl_factory = gtk_icon_factory_new();

	respath = g_strdup_printf("%sgretl-icons.gresource", gretl_home());
	icons = g_resource_load(respath, NULL);
	if (icons == NULL) {
	    fprintf(stderr, "g_resource_load: failed to load icons\n");
	    g_free(respath);
	    goto do_pixmaps;
	}

	g_resources_register(icons);

	if (bigger) {
	    strcpy(icon_path, "/gretl/icons/24x24/");
	    strcpy(menu_path, "/gretl/icons/16x16/");
	    pm = strrchr(menu_path, '/') + 1;
	} else {
	    strcpy(icon_path, "/gretl/icons/16x16/");
	}
	p = strrchr(icon_path, '/') + 1;

	for (i=0; i<n1; i++) {
	    strcat(icon_path, png_stocks[i].fname);
	    pbuf = gdk_pixbuf_new_from_resource(icon_path, NULL);
	    if (pbuf == NULL) {
		fprintf(stderr, "Failed to load %s\n", icon_path);
		*p = '\0';
		continue;
	    }
	    if (bigger && png_stocks[i].in_menu) {
		iset = gtk_icon_set_new();
		/* for toolbar use */
		isrc = gtk_icon_source_new();
		gtk_icon_source_set_pixbuf(isrc, pbuf);
		gtk_icon_source_set_size(isrc, toolbar_icon_size);
		gtk_icon_source_set_size_wildcarded(isrc, FALSE);
		gtk_icon_set_add_source(iset, isrc);
		g_object_unref(pbuf);
		/* for menu use */
		strcat(menu_path, png_stocks[i].fname);
		pbuf = gdk_pixbuf_new_from_resource(menu_path, NULL);
		isrc = gtk_icon_source_new();
		gtk_icon_source_set_pixbuf(isrc, pbuf);
		gtk_icon_source_set_size(isrc, GTK_ICON_SIZE_MENU);
		gtk_icon_source_set_size_wildcarded(isrc, FALSE);
		gtk_icon_set_add_source(iset, isrc);
		g_object_unref(pbuf);
		*pm = '\0';
	    } else {
		/* we just need a single icon */
		iset = gtk_icon_set_new_from_pixbuf(pbuf);
	    }
	    gtk_icon_factory_add(gretl_factory, png_stocks[i].id, iset);
	    gtk_icon_set_unref(iset);
	    *p = '\0';
	}

	g_free(respath);
	g_resources_unregister(icons);
	g_resource_unref(icons);

    do_pixmaps:

	for (i=0; i<n2; i++) {
	    pbuf = gdk_pixbuf_new_from_xpm_data((const char **) xpm_stocks[i].xpm);
	    iset = gtk_icon_set_new_from_pixbuf(pbuf);
	    g_object_unref(pbuf);
	    gtk_icon_factory_add(gretl_factory, xpm_stocks[i].id, iset);
	    gtk_icon_set_unref(iset);
	}

	gtk_icon_factory_add_default(gretl_factory);
    }
}

static void save_as_callback (GtkWidget *w, windata_t *vwin)
{
    GtkWidget *vmain = vwin_toplevel(vwin);
    int ci = 0;

    if (g_object_get_data(G_OBJECT(vmain), "text_out")) {
	ci = SAVE_OUTPUT;
    } else if (vwin->role == EDIT_HANSL) {
	ci = SAVE_SCRIPT;
    } else if (vwin->role == EDIT_GP) {
	ci = SAVE_GP_CMDS;
    } else if (vwin->role == EDIT_R) {
	ci = SAVE_R_CMDS;
    } else if (vwin->role == EDIT_OX) {
	ci = SAVE_OX_CMDS;
    } else if (vwin->role == EDIT_OCTAVE) {
	ci = SAVE_OCTAVE_CMDS;
    } else if (vwin->role == EDIT_PYTHON) {
        ci = SAVE_PYTHON_CMDS;
    } else if (vwin->role == EDIT_JULIA) {
	ci = SAVE_JULIA_CODE;
    } else if (vwin->role == EDIT_DYNARE) {
	ci = SAVE_DYNARE_CODE;
    } else if (vwin->role == EDIT_LPSOLVE) {
	ci = SAVE_LPSOLVE_CODE;
    } else if (vwin->role == EDIT_STATA) {
	ci = SAVE_STATA_CMDS;
    } else if (vwin->role == EDIT_SPEC) {
	ci = SAVE_SPEC_FILE;
    } else if (vwin->role == VIEW_FILE) {
	ci = SAVE_TEXT;
    } else if (vwin->role == EDIT_X12A) {
	ci = SAVE_X13_SPC;
    } else {
	dummy_call();
	return;
    }

    if (ci == SAVE_TEXT) {
	file_selector(ci, FSEL_DATA_MISC, vwin->data);
    } else {
	file_selector(ci, FSEL_DATA_VWIN, vwin);
    }
}

/* callback for the "Open" icon in a script editing window,
   which enables the user to switch to a different script,
   or to open another tab if the editor is tab-enabled
*/

static void file_open_callback (GtkWidget *w, windata_t *vwin)
{
    file_selector(OPEN_SCRIPT, FSEL_DATA_VWIN, vwin);
}

static void toolbar_new_callback (GtkWidget *w, windata_t *vwin)
{
    do_new_script(vwin->role, NULL, NULL);
}

static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    if (textview_use_highlighting(vwin->role)) {
	int resp = yes_no_cancel_dialog(NULL,
					_("Print with syntax highlighting?"),
					vwin_toplevel(vwin));

	if (resp == GRETL_YES) {
	    sourceview_print(vwin);
	} else if (resp == GRETL_NO) {
	    window_print(NULL, vwin);
	}
    } else {
	window_print(NULL, vwin);
    }
}

static void split_pane_callback (GtkWidget *w, windata_t *vwin)
{
    GtkWidget *hb = g_object_get_data(G_OBJECT(w), "hpane");
    GtkWidget *vb = g_object_get_data(G_OBJECT(w), "vpane");
    int vertical = 0;

    if (hb != NULL) {
	vb = w;
	vertical = 1;
    } else {
	hb = w;
    }

    /* Note: by "vertical" here we mean that the split runs vertically,
       dividing the pane into left- and right-hand sections; otherwise
       the split runs horizontally. In a "gnuplot commands" window we
       only offer a horizontal split, which means that @vb ("vertical
       button") may be NULL.
    */

    if (g_object_get_data(G_OBJECT(vwin->vbox), "sw") != NULL) {
	/* currently in single-view mode: so split */
	viewer_split_pane(vwin, vertical);
	if (vb != NULL) {
	    gtk_widget_set_sensitive(vb, vertical);
	}
	gtk_widget_set_sensitive(hb, !vertical);
	if (vertical) {
	    gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(vb),
					 GRETL_STOCK_JOIN_V);
	} else {
	    gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(hb),
					 GRETL_STOCK_JOIN_H);
	}
    } else {
	GtkWidget *paned;

	paned = g_object_get_data(G_OBJECT(vwin->vbox), "paned");

	if (paned != NULL) {
	    /* currently in split-view mode: so rejoin */
	    vertical = GTK_IS_HPANED(paned);
	    viewer_close_pane(vwin);
	    gtk_widget_set_sensitive(hb, TRUE);
	    if (vb != NULL) {
		gtk_widget_set_sensitive(vb, TRUE);
	    }
	    if (vertical) {
		gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(vb),
					     GRETL_STOCK_SPLIT_V);
	    } else {
		gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(hb),
					     GRETL_STOCK_SPLIT_H);
	    }
	}
    }
}

static void editor_prefs_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == CONSOLE) {
	console_prefs_dialog(vwin->main);
    } else {
	preferences_dialog(TAB_EDITOR, NULL, vwin_toplevel(vwin));
    }
}

static void activate_script_help (GtkWidget *widget, windata_t *vwin)
{
    text_set_cursor(vwin->text, GDK_QUESTION_ARROW);
    set_window_help_active(vwin);
}

static int edit_script_popup_item (GretlToolItem *item)
{
    if (item->icon == NULL) return 0;

    return !strcmp(item->icon, GTK_STOCK_COPY) ||
	!strcmp(item->icon, GTK_STOCK_PASTE) ||
	!strcmp(item->icon, GTK_STOCK_FIND) ||
	!strcmp(item->icon, GTK_STOCK_UNDO) ||
	!strcmp(item->icon, GTK_STOCK_FIND_AND_REPLACE);
}

static void vwin_cut_callback (GtkWidget *w, windata_t *vwin)
{
    gtk_text_buffer_cut_clipboard(gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text)),
				  gtk_clipboard_get(GDK_NONE),
				  TRUE);
}

static GretlToolItem viewbar_items[] = {
    { N_("New window"), GTK_STOCK_NEW, G_CALLBACK(toolbar_new_callback), NEW_ITEM },
    { N_("Open..."), GTK_STOCK_OPEN, G_CALLBACK(file_open_callback), OPEN_ITEM },
    { N_("Save"), GTK_STOCK_SAVE, G_CALLBACK(vwin_save_callback), SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, G_CALLBACK(save_as_callback), SAVE_AS_ITEM },
    { N_("Print..."), GTK_STOCK_PRINT, G_CALLBACK(window_print_callback), PRINT_ITEM },
    { N_("Run"), GTK_STOCK_EXECUTE, G_CALLBACK(do_run_script), EXEC_ITEM },
    { N_("Cut"), GTK_STOCK_CUT, G_CALLBACK(vwin_cut_callback), EDIT_ITEM },
    { N_("Copy"), GTK_STOCK_COPY, G_CALLBACK(vwin_copy_callback), COPY_ITEM },
    { N_("Paste"), GTK_STOCK_PASTE, G_CALLBACK(text_paste), EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, G_CALLBACK(text_find), FIND_ITEM },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, G_CALLBACK(text_replace), EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, G_CALLBACK(text_undo), EDIT_ITEM },
    { N_("Redo"), GTK_STOCK_REDO, G_CALLBACK(text_redo), EDIT_ITEM },
    { N_("Preferences..."), GTK_STOCK_PREFERENCES, G_CALLBACK(editor_prefs_callback), EDIT_HANSL_ITEM },
    { N_("Auto-indent script"), GTK_STOCK_INDENT, G_CALLBACK(indent_hansl), EDIT_HANSL_ITEM },
    { N_("Toggle split pane"), GRETL_STOCK_SPLIT_H, G_CALLBACK(split_pane_callback), SPLIT_H_ITEM },
    { N_("Toggle split pane"), GRETL_STOCK_SPLIT_V, G_CALLBACK(split_pane_callback), SPLIT_V_ITEM },
    { N_("Help on command"), GRETL_STOCK_QUERY, G_CALLBACK(activate_script_help), CMD_HELP_ITEM }
};

static int n_viewbar_items = G_N_ELEMENTS(viewbar_items);

#define exec_ok(r) (vwin_editing_script(r))
#define open_ok(r) (vwin_editing_script(r))
#define new_ok(r)  (vwin_editing_script(r))
#define edit_ok(r) (vwin_editing_script(r))

#define save_as_ok(r) (r != EDIT_HEADER && \
	               r != EDIT_NOTES && \
	               r != EDIT_PKG_CODE && \
		       r != EDIT_PKG_SAMPLE && \
		       r != CONSOLE && \
		       r != VIEW_BUNDLE && \
		       r != VIEW_DBNOMICS)

#define cmd_help_ok(r) (r == EDIT_HANSL)

#define split_h_ok(r) (r == SCRIPT_OUT || vwin_editing_script(r))

#define split_v_ok(r) (r == SCRIPT_OUT)

/* Screen out unwanted menu items depending on the context; also
   adjust the callbacks associated with some items based on
   context.
*/

static GCallback tool_item_get_callback (GretlToolItem *item, windata_t *vwin,
					 int save_ok)
{
    GCallback func = item->func;
    int f = item->flag;
    int r = vwin->role;

    if (use_toolbar_search_box(r) && f == FIND_ITEM) {
	/* using an "inline" search box: skip the
	   "Find" button */
	return NULL;
    }

    if (!edit_ok(r) && f == EDIT_ITEM) {
	return NULL;
    } else if (!open_ok(r) && f == OPEN_ITEM) {
	return NULL;
    } else if (!new_ok(r) && f == NEW_ITEM) {
	return NULL;
    } else if (!exec_ok(r) && f == EXEC_ITEM) {
	return NULL;
    } else if (!cmd_help_ok(r) && f == CMD_HELP_ITEM) {
	return NULL;
    } else if (!split_h_ok(r) && f == SPLIT_H_ITEM) {
	return NULL;
    } else if (!split_v_ok(r) && f == SPLIT_V_ITEM) {
	return NULL;
    } else if (r != EDIT_HANSL && f == EDIT_HANSL_ITEM) {
	return NULL;
    } else if (f == SAVE_ITEM && !save_ok) {
	return NULL;
    } else if (f == SAVE_AS_ITEM) {
	if (!save_as_ok(r) || (vwin->flags & VWIN_NO_SAVE)) {
	    return NULL;
	}
    }

    return func;
}

static void gretl_toolbar_flat (GtkWidget *w)
{
    static int style_done;

    gtk_widget_set_name(w, "gretl_toolbar");

    if (!style_done) {
	gtk_rc_parse_string("style \"gretl-tb-style\"\n{\n"
			    "  GtkToolbar::shadow-type = GTK_SHADOW_NONE\n"
			    "}\n"
			    "widget \"*.gretl_toolbar\" style \"gretl-tb-style\"");
	style_done = 1;
    }
}

GtkWidget *gretl_toolbar_new (GtkWidget *sibling)
{
    GtkWidget *tb = gtk_toolbar_new();

    gtk_toolbar_set_icon_size(GTK_TOOLBAR(tb), toolbar_icon_size);
    gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_show_arrow(GTK_TOOLBAR(tb), FALSE);

    if (sibling == NULL) {
	/* if we're not alongside a menu bar ("sibling"),
	   show the toolbar without a shadow
	*/
	gretl_toolbar_flat(tb);
    }

    return tb;
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_widget_set_tooltip_text(w, str);
}

static GtkToolItem *gretl_menu_button (const char *icon,
				       const char *tip,
				       GtkWidget **pw)
{
    GtkWidget *img, *button = gtk_button_new();
    GtkToolItem *item = gtk_tool_item_new();

    gtk_widget_set_tooltip_text(GTK_WIDGET(item), _(tip));
    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    img = gtk_image_new_from_stock(icon, toolbar_icon_size /* GTK_ICON_SIZE_MENU */);
    gtk_container_add(GTK_CONTAINER(button), img);
    gtk_container_add(GTK_CONTAINER(item), button);
    *pw = button;

    return item;
}

static void gretl_tool_item_set_tip (GtkWidget *item,
				     GretlToolItem *tool)
{
    const char *accel = NULL;

    if (tool->flag == EXEC_ITEM) {
	accel = "Ctrl+R";
    } else if (tool->flag == COPY_ITEM) {
	accel = "Ctrl+C";
    } else if (tool->flag == SAVE_ITEM) {
	accel = "Ctrl+S";
    } else if (tool->flag == FIND_ITEM) {
	accel = "Ctrl+F";
    } else if (!strcmp(tool->icon, GTK_STOCK_FIND_AND_REPLACE)) {
	accel = "Ctrl+H";
    }

    if (accel != NULL) {
	gchar *s = g_strdup_printf("%s (%s)", _(tool->tip), accel);

	gtk_widget_set_tooltip_text(item, s);
	g_free(s);
    } else {
	gtk_widget_set_tooltip_text(item, _(tool->tip));
    }
}

GtkWidget *gretl_toolbar_insert (GtkWidget *tbar,
				 GretlToolItem *tool,
				 GCallback func,
				 gpointer data,
				 gint pos)
{
    GtkToolItem *item;

    item = gtk_tool_button_new_from_stock(tool->icon);
    gretl_tool_item_set_tip(GTK_WIDGET(item), tool);
    g_signal_connect(G_OBJECT(item), "clicked", func, data);
    gtk_widget_set_size_request(GTK_WIDGET(item), 30, -1);
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), item, pos);

    return GTK_WIDGET(item);
}

static void button_menu_pos (GtkMenu *menu,
			     gint *x,
			     gint *y,
			     gboolean *push_in,
			     gpointer data)
{
    GtkWidget *button = data;
    gint wx, wy, tx, ty;

    gdk_window_get_origin(gtk_widget_get_window(button), &wx, &wy);
    gtk_widget_translate_coordinates(button, gtk_widget_get_toplevel(button),
				     0, 0, &tx, &ty);
    *x = wx + tx;
    *y = wy + ty + 26;
    *push_in = TRUE;
}

static void tool_item_popup (GtkWidget *button, GdkEvent *event,
			     GtkWidget *menu)
{
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL,
		   button_menu_pos, button,
		   event->button.button, event->button.time);
}

GtkWidget *vwin_toolbar_insert (GretlToolItem *tool,
				GCallback func,
				GtkWidget *menu,
				windata_t *vwin,
				gint pos)
{
    GtkToolItem *item;

    if (menu != NULL) {
	/* make and insert a button that pops down a menu */
	GtkWidget *button;

	item = gretl_menu_button(tool->icon, tool->tip, &button);
	g_signal_connect(G_OBJECT(button), "button-press-event",
			 G_CALLBACK(tool_item_popup), menu);
    } else {
	/* make and insert a regular callback button */
	item = gtk_tool_button_new_from_stock(tool->icon);
	g_signal_connect(G_OBJECT(item), "clicked", func, vwin);
	if (tool->flag == NEW_ITEM && window_is_tab(vwin)) {
	    gtk_widget_set_tooltip_text(GTK_WIDGET(item), _("New tab"));
	} else {
	    gretl_tool_item_set_tip(GTK_WIDGET(item), tool);
	}
    }

    gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), item, pos);

    return GTK_WIDGET(item);
}

static void viewbar_add_items (windata_t *vwin, ViewbarFlags flags)
{
    int save_ok = (flags & VIEWBAR_EDITABLE);
    GtkWidget *hpane = NULL, *vpane = NULL;
    GtkWidget *button;
    GtkWidget *menu;
    GretlToolItem *item;
    GCallback func;
    int i;

    for (i=0; i<n_viewbar_items; i++) {
	func = NULL;
	menu = NULL;
	item = &viewbar_items[i];

	/* Is there anything to hook up, in context? We
	   try first for a menu to attach to the toolbar
	   button; failing that we test for a "direct"
	   callback function.
	*/
	//menu = tool_item_get_menu(item, vwin);
	if (menu == NULL && item->func != NULL) {
	    func = tool_item_get_callback(item, vwin, save_ok);
	}
	if (func == NULL && menu == NULL) {
	    /* nothing to hook up */
	    continue;
	}

	button = vwin_toolbar_insert(item, func, menu, vwin, -1);

	if (func == (GCallback) split_pane_callback) {
	    if (hpane == NULL) {
		hpane = button;
	    } else {
		vpane = button;
	    }
	}

	if (item->flag == SAVE_ITEM) {
	    /* nothing to save just yet */
	    g_object_set_data(G_OBJECT(vwin->mbar), "save_button", button);
	    gtk_widget_set_sensitive(button, FALSE);
	} else if (item->flag == SAVE_AS_ITEM) {
	    g_object_set_data(G_OBJECT(vwin->mbar), "save_as_button", button);
	    if (strstr(vwin->fname, "script_tmp")) {
		gtk_widget_set_sensitive(button, FALSE);
	    }
	}
    }

    if (hpane != NULL) {
	g_object_set_data(G_OBJECT(hpane), "vpane", vpane);
    }
    if (vpane != NULL) {
	g_object_set_data(G_OBJECT(vpane), "hpane", hpane);
    }
}

void vwin_add_viewbar (windata_t *vwin, ViewbarFlags flags)
{
    if ((flags & VIEWBAR_HAS_TEXT) || vwin->role == SCRIPT_OUT) {
	g_object_set_data(G_OBJECT(vwin->main), "text_out",
			  GINT_TO_POINTER(1));
    }

    vwin->mbar = gretl_toolbar_new(NULL);
    viewbar_add_items(vwin, flags);
    vwin_pack_toolbar(vwin);
}

GtkWidget *build_text_popup (windata_t *vwin)
{
    GtkWidget *pmenu = gtk_menu_new();
    GretlToolItem *item;
    GCallback func;
    GtkWidget *w;
    int i;

    for (i=0; i<n_viewbar_items; i++) {
	item = &viewbar_items[i];
	func = G_CALLBACK(NULL);
	if (item->flag == SPLIT_H_ITEM || item->flag == SPLIT_V_ITEM) {
	    continue;
	} else if (vwin->role == EDIT_HANSL) {
	    /* the script editor popup may have some special stuff
	       added: don't clutter it up */
	    if (edit_script_popup_item(item)) {
		func = item->func;
	    } else {
		func = NULL;
	    }
	} else {
	    func = tool_item_get_callback(item, vwin, 0);
	}
	if (func != G_CALLBACK(NULL)) {
	    if (func == G_CALLBACK(text_paste)) {
		GtkClipboard *cb = gtk_clipboard_get(GDK_NONE);

		if (!gtk_clipboard_wait_is_text_available(cb)) {
		    continue;
		}
	    } else if (func == G_CALLBACK(text_undo) && !text_can_undo(vwin)) {
		continue;
	    }
	    w = gtk_menu_item_new_with_label(_(item->tip));
	    g_signal_connect(G_OBJECT(w), "activate", func, vwin);
	    gtk_widget_show(w);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), w);
	}
    }

    return pmenu;
}

/* Add a temporary menubar for use in a script output
   window, while we're waiting for the output. If the
   output window is being reused this is a bit more
   complicated; we have to "hide" the regular menubar
   before inserting the temporary one.
 */

void vwin_add_tmpbar (windata_t *vwin)
{
    GretlToolItem item = {
	N_("Stop"),
	GTK_STOCK_STOP,
	G_CALLBACK(do_stop_script),
	0
    };
    GtkWidget *hbox, *tmp;

    hbox = g_object_get_data(G_OBJECT(vwin->main), "top-hbox");

    if (hbox != NULL) {
	/* We're replacing a "real" menubar temporarily: ref. the
	   widgets in @hbox before removing them so we can put
	   them back later.
	*/
	GtkWidget *winlist = g_object_get_data(G_OBJECT(hbox), "winlist");

	g_object_ref(G_OBJECT(vwin->mbar));
	gtk_container_remove(GTK_CONTAINER(hbox), vwin->mbar);
	if (vwin->finder != NULL) {
	    g_object_ref(G_OBJECT(vwin->finder));
	    gtk_container_remove(GTK_CONTAINER(hbox), vwin->finder);
	}
	if (winlist != NULL) {
	    g_object_ref(G_OBJECT(winlist));
	    gtk_container_remove(GTK_CONTAINER(hbox), winlist);
	}
    } else {
	/* starting from scratch */
	hbox = gtk_hbox_new(FALSE, 0);
	g_object_set_data(G_OBJECT(vwin->main), "top-hbox", hbox);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
    }

    tmp = gretl_toolbar_new(NULL);
    gretl_toolbar_insert(tmp, &item, item.func, NULL, 0);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    start_wait_for_output(vwin, hbox);
    gtk_widget_show_all(hbox);
}
