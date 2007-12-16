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

#define NEWSTYLE 1

#if (NEWSTYLE == 0)
/* pixmaps for gretl toolbar */
#include "../pixmaps/mini.calc.xpm"
#include "../pixmaps/mini.edit.xpm"
#include "../pixmaps/mini.sh.xpm"
#include "../pixmaps/mini.session.xpm"
#include "../pixmaps/mini.manual.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.plot.xpm"
#include "../pixmaps/mini.model.xpm"
#include "../pixmaps/mini.ofolder.xpm"
#include "../pixmaps/mini.browser.xpm"
#endif

static GtkWidget *toolbar_box;

/* callbacks for gretl toolbar icons */

static void show_calc (void)
{
#ifdef G_OS_WIN32
    create_child_process(calculator);
#else
    gretl_fork(calculator, NULL);
#endif 
}

static void open_textbook_data (void)
{
    display_files(NULL, TEXTBOOK_DATA, NULL);
}

static void gretl_website (void)
{
    if (browser_open("http://gretl.sourceforge.net/")) {
	errbox("Failed to open URL");
    }
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
    do_new_script(NULL, 0, NULL);
}

/* end toolbar icon callbacks */

#if NEWSTYLE

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
    { N_("gretl website"), GRETL_STOCK_WWW, gretl_website },
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

    toolbar_box = gtk_handle_box_new();
    gtk_box_pack_start(GTK_BOX(hbox), toolbar_box, FALSE, FALSE, 0);

    toolbar = gtk_toolbar_new();
    gtk_container_add(GTK_CONTAINER(toolbar_box), toolbar);

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

#else

static GtkWidget *image_button_new (GdkPixbuf *pix, void (*toolfunc)())
{
    GtkWidget *image = gtk_image_new_from_pixbuf(pix);
    GtkWidget *button = gtk_button_new();

    gtk_widget_set_size_request(button, 26, 24);

    gtk_container_add(GTK_CONTAINER(button), image);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(toolfunc), NULL);

    return button;
}

static void make_toolbar (GtkWidget *vbox)
{
    GtkWidget *button;
    GtkWidget *toolbar;
    GtkWidget *hbox;
    GdkPixbuf *icon;
    int i;
    const char *toolstrings[] = {
	N_("launch calculator"), 
	N_("new script"), 
	N_("open gretl console"),
	N_("session icon view"),
	N_("gretl website"), 
	N_("user's guide"),
	N_("command reference"), 
	N_("X-Y graph"), 
	N_("OLS model"),
	N_("open dataset"),
	NULL
    };
    gchar **toolxpm = NULL;
    void (*toolfunc)() = NULL;
    const char *toolstr;

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    toolbar_box = gtk_handle_box_new();
    gtk_box_pack_start(GTK_BOX(hbox), toolbar_box, FALSE, FALSE, 0);

    toolbar = gtk_hbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(toolbar_box), toolbar);

    for (i=0; toolstrings[i] != NULL; i++) {
	switch (i) {
	case 0:
	    toolxpm = mini_calc_xpm;
	    toolfunc = show_calc;
	    break;
	case 1:
	    toolxpm = mini_edit_xpm;
	    toolfunc = new_script_callback;
	    break;
	case 2:
	    toolxpm = mini_sh_xpm;
	    toolfunc = show_gretl_console;
	    break;
	case 3:
	    toolxpm = mini_session_xpm;
	    toolfunc = go_session;
	    break;
	case 4:
	    toolxpm = mini_browser_xpm;
	    toolfunc = gretl_website;
	    break;  
	case 5:
	    toolxpm = mini_pdf_xpm;
	    toolfunc = toolbar_users_guide;
	    break;    
	case 6:
	    toolxpm = mini_manual_xpm;
	    toolfunc = toolbar_command_reference;
	    break;
	case 7:
	    toolxpm = mini_plot_xpm;
	    toolfunc = xy_graph;
	    break;
	case 8:
	    toolxpm = mini_model_xpm;
	    toolfunc = ols_model;
	    break;
	case 9:
	    toolxpm = mini_ofolder_xpm;
	    toolfunc = open_textbook_data;
	    break;
	default:
	    break;
	}

	toolstr = _(toolstrings[i]);
	icon = gdk_pixbuf_new_from_xpm_data((const char **) toolxpm);
	button = image_button_new(icon, toolfunc);
	gtk_box_pack_start(GTK_BOX(toolbar), button, FALSE, FALSE, 0);
	gretl_tooltips_add(button, toolstr);
	gdk_pixbuf_unref(icon);
    }

    gtk_widget_show_all(hbox);
}

#endif

/* public interface */

void show_or_hide_toolbar (int want_toolbar)
{
    if (want_toolbar && toolbar_box == NULL) {
	GtkWidget *vbox = g_object_get_data(G_OBJECT(mdata->w), "vbox");

	make_toolbar(vbox);
    } else if (!want_toolbar && toolbar_box != NULL) {
	gtk_widget_destroy(toolbar_box);
	toolbar_box = NULL;
    }
}

