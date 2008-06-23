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

/* toolbar.c: the gretl toolbar */

#include "gretl.h"
#include "console.h"
#include "session.h"
#include "datafiles.h"
#include "selector.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

/* callbacks for gretl toolbar icons */

static void tbar_calc (void)
{
#ifdef G_OS_WIN32
    create_child_process(calculator);
#else
    gretl_fork("calculator", NULL);
#endif 
}

static void tbar_open_data (void)
{
    display_files(TEXTBOOK_DATA, NULL);
}

static void tbar_users_guide (void)
{
    display_pdf_help(NULL);
}

static void tbar_command_ref (void)
{
    plain_text_cmdref(NULL);
}

static void tbar_xy_graph (void)
{
    if (data_status) {
	if (datainfo->v == 2) {
	    do_graph_var(mdata->active_var);
	} else if (mdata_selection_count() == 2) {
	    plot_from_selection(GR_XY);
	} else {
	    selection_dialog(_("gretl: define graph"), 
			     do_graph_from_selector,
			     GR_XY, 0);
	}
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_model (void)
{
    if (data_status) {
	selection_dialog(_("gretl: specify model"), do_model, OLS, 0);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_iconview (void)
{
    if (data_status) {
	view_session();
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_new_script (void)
{
    do_new_script(EDIT_SCRIPT);
}

static void tbar_show_funcs (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, w);
}

/* end toolbar icon callbacks */

struct toolbar_item {
    const char *tip;
    const gchar *icon;
    void (*toolfunc)();
};

#if NO_EDIT_ICON
# define GTK_STOCK_EDIT "gretl-script"
#endif

static struct toolbar_item toolbar_items[] = {
    { N_("launch calculator"),  GRETL_STOCK_CALC,    tbar_calc },
    { N_("new script"),         GTK_STOCK_EDIT,      tbar_new_script },
    { N_("open gretl console"), GRETL_STOCK_CONSOLE, show_gretl_console },
    { N_("session icon view"),  GRETL_STOCK_ICONS,   tbar_iconview },
    { N_("function packages"),  GRETL_STOCK_FUNC,    tbar_show_funcs },
    { N_("user's guide"),       GRETL_STOCK_PDF,     tbar_users_guide },
    { N_("command reference"),  GTK_STOCK_HELP,      tbar_command_ref },
    { N_("X-Y graph"),          GRETL_STOCK_SCATTER, tbar_xy_graph },
    { N_("OLS model"),          GRETL_STOCK_MODEL,   tbar_model },
    { N_("open dataset"),       GTK_STOCK_OPEN,      tbar_open_data }
};

static void make_toolbar (GtkWidget *vbox)
{
    GtkWidget *hbox, *toolbar;
    GtkToolItem *item;
    int i, n = G_N_ELEMENTS(toolbar_items);

    gretl_stock_icons_init();

    toolbar = gtk_toolbar_new();
    gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar), GTK_ICON_SIZE_MENU);
    gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_show_arrow(GTK_TOOLBAR(toolbar), FALSE);

    for (i=0; i<n; i++) {
	item = gtk_tool_button_new_from_stock(toolbar_items[i].icon);
	gtk_tool_item_set_tooltip(item, get_gretl_tips(), 
				  _(toolbar_items[i].tip),
				  NULL);
	g_signal_connect(item, "clicked", toolbar_items[i].toolfunc, mdata);
	gtk_toolbar_insert(GTK_TOOLBAR(toolbar), item, -1);
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), toolbar, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

/* public interface */

void show_toolbar (void)
{
    GtkWidget *vbox = g_object_get_data(G_OBJECT(mdata->w), "vbox");

    make_toolbar(vbox);
}

