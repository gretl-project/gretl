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
    display_files(NULL, TEXTBOOK_DATA, NULL);
}

void toolbar_users_guide (void)
{
    display_pdf_help(NULL, 1, NULL);
}

void toolbar_command_reference (void)
{
    plain_text_cmdref(NULL, 0, NULL);
}

static void xy_graph (void)
{
    if (data_status) {
	if (datainfo->v == 2) {
	    do_graph_var(mdata->active_var);
	} else if (mdata_selection_count() == 2) {
	    plot_from_selection(NULL, GR_XY, NULL);
	} else {
	    selector_callback(NULL, GR_XY, NULL);
	}
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void ols_model (void)
{
    if (data_status) {
	model_callback(NULL, OLS, NULL);
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

static void new_script_callback (void)
{
    do_new_script(NULL, EDIT_SCRIPT, NULL);
}

static void show_funcs_callback (GtkWidget *w, gpointer p)
{
    display_files(NULL, FUNC_FILES, w);
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
    { N_("new script"), GRETL_STOCK_SCRIPT, new_script_callback },
#else
    { N_("new script"), GTK_STOCK_EDIT, new_script_callback },
#endif
    { N_("open gretl console"), GRETL_STOCK_CONSOLE, show_gretl_console },
    { N_("session icon view"), GRETL_STOCK_ICONS, go_session },
    { N_("function packages"), GRETL_STOCK_FUNC, show_funcs_callback },
    { N_("user's guide"), GRETL_STOCK_PDF, toolbar_users_guide },
    { N_("command reference"), GTK_STOCK_HELP, toolbar_command_reference },
    { N_("X-Y graph"), GRETL_STOCK_SCATTER, xy_graph },
    { N_("OLS model"), GRETL_STOCK_MODEL, ols_model },
    { N_("open dataset"), GTK_STOCK_OPEN, open_textbook_data },
    { NULL, NULL, NULL }
};

static void make_toolbar (GtkWidget *vbox)
{
    GtkWidget *w, *hbox;
    GtkWidget *toolbar, *image;
    int i;

    gretl_stock_icons_init();

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    toolbar = gtk_toolbar_new();
    gtk_box_pack_start(GTK_BOX(hbox), toolbar, FALSE, FALSE, 0);

    for (i=0; toolbar_items[i].str != NULL; i++) {
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

