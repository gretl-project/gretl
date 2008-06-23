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

static void show_calc (void)
{
#ifdef G_OS_WIN32
    create_child_process(calculator);
#else
    gretl_fork("calculator", NULL);
#endif 
}

static void open_textbook_data (void)
{
    display_files(TEXTBOOK_DATA, NULL);
}

static void toolbar_users_guide (void)
{
    display_pdf_help(NULL);
}

static void toolbar_command_reference (void)
{
    plain_text_cmdref(NULL);
}

static void xy_graph (void)
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

static void ols_model (void)
{
    if (data_status) {
	selection_dialog(_("gretl: specify model"), do_model, OLS, 0);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void go_session (void)
{
    if (data_status) {
	view_session();
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void toolbar_new_script (void)
{
    do_new_script(EDIT_SCRIPT);
}

static void show_funcs_callback (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, w);
}

/* end toolbar icon callbacks */

struct toolbar_item {
    const char *str;
    const gchar *icon;
    void (*toolfunc)();
};

static struct toolbar_item toolbar_items[] = {
    { N_("launch calculator"), GRETL_STOCK_CALC, show_calc },
#if NO_EDIT_ICON
    { N_("new script"), GRETL_STOCK_SCRIPT, toolbar_new_script },
#else
    { N_("new script"), GTK_STOCK_EDIT, toolbar_new_script },
#endif
    { N_("open gretl console"), GRETL_STOCK_CONSOLE, show_gretl_console },
    { N_("session icon view"), GRETL_STOCK_ICONS, go_session },
    { N_("function packages"), GRETL_STOCK_FUNC, show_funcs_callback },
    { N_("user's guide"), GRETL_STOCK_PDF, toolbar_users_guide },
    { N_("command reference"), GTK_STOCK_HELP, toolbar_command_reference },
    { N_("X-Y graph"), GRETL_STOCK_SCATTER, xy_graph },
    { N_("OLS model"), GRETL_STOCK_MODEL, ols_model },
    { N_("open dataset"), GTK_STOCK_OPEN, open_textbook_data }
};

#if 0 /* not yet */
static GtkActionEntry toolbar_items[] = {
    { "calc", GRETL_STOCK_CALC, N_("launch calculator"), show_calc },
#if NO_EDIT_ICON
    { "newscript", GRETL_STOCK_SCRIPT, N_("new script"), toolbar_new_script },
#else
    { "newscript", GTK_STOCK_EDIT, N_("new script"), toolbar_new_script },
#endif
    { "console", GRETL_STOCK_CONSOLE, N_("open gretl console"), 
      show_gretl_console },
    { "iconview", GRETL_STOCK_ICONS, N_("session icon view"), go_session },
    { "fnpkgs", GRETL_STOCK_FUNC, N_("function packages"), show_funcs_callback },
    { "userguide", GRETL_STOCK_PDF, N_("user's guide"), toolbar_users_guide },
    { "cmdref", GTK_STOCK_HELP, N_("command reference"), 
      toolbar_command_reference },
    { "xygraph", GRETL_STOCK_SCATTER, N_("X-Y graph"), xy_graph },
    { "ols", GRETL_STOCK_MODEL, N_("OLS model"), ols_model },
    { "opendata", GTK_STOCK_OPEN, N_("open dataset"), open_textbook_data }
};
#endif

static void make_toolbar (GtkWidget *vbox)
{
    GtkWidget *w, *hbox;
    GtkWidget *toolbar, *image;
    int i, n = G_N_ELEMENTS(toolbar_items);

    gretl_stock_icons_init();

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    toolbar = gtk_toolbar_new();
    gtk_box_pack_start(GTK_BOX(hbox), toolbar, FALSE, FALSE, 0);

    for (i=0; i<n; i++) {
	image = gtk_image_new();
	gtk_image_set_from_stock(GTK_IMAGE(image), toolbar_items[i].icon, 
				 GTK_ICON_SIZE_MENU);
        w = gtk_toolbar_append_item(GTK_TOOLBAR(toolbar),
				    NULL, _(toolbar_items[i].str), NULL,
				    image, toolbar_items[i].toolfunc, mdata);
    }

    gtk_widget_show_all(hbox);
}

/* public interface */

void show_toolbar (void)
{
    GtkWidget *vbox = g_object_get_data(G_OBJECT(mdata->w), "vbox");

    make_toolbar(vbox);
}

