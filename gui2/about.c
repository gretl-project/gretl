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
#include "version.h"
#include "dlgutils.h"
#include "build.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

const gchar *copyright = "Copyright (C) 2000-2014 Allin Cottrell and "
                         "Riccardo \"Jack\" Lucchetti";
const gchar *bonmot = N_("\"By econometricians, for econometricians.\"");
const gchar *website = "http://gretl.sourceforge.net/";

static GtkWidget *open_logo (void)
{
    char *fname;
    GdkPixbuf *pbuf;
    GError *error = NULL;
    GtkWidget *image = NULL;

    fname = g_strdup_printf("%sgretl-logo.xpm", gretl_home());
    pbuf = gdk_pixbuf_new_from_file(fname, &error);

    if (pbuf == NULL) {
	errbox(error->message);
	g_error_free(error);
    } else {
	image = gtk_image_new_from_pixbuf(pbuf);
	g_object_unref(pbuf);
    }

    g_free(fname);

    return image;
}

static void license_callback (GtkWidget *w, gpointer p)
{
    gchar *fname;

    fname = g_strdup_printf("%sCOPYING", gretl_home());
    view_file(fname, 0, 0, 78, 350, VIEW_FILE);
    g_free(fname);
}

static void relnotes_callback (GtkWidget *w, gpointer p)
{
    gchar *fname;

    fname = g_strdup_printf("%sNEWS", gretl_home());
    view_file(fname, 0, 0, 78, 350, VIEW_FILE);
    g_free(fname);
}

static void show_link_cursor (GtkWidget *w, gpointer p)
{
    GdkWindow *window = gtk_widget_get_window(w);
    GdkCursor *c;

    c = gdk_cursor_new(GDK_HAND2);
    gdk_window_set_cursor(window, c);
    gdk_cursor_unref(c);
}

static void show_website (GtkWidget *w, gpointer p)
{
    if (browser_open(website)) {
	errbox("Failed to open URL");
    }
}

void about_dialog (void) 
{
    GtkWidget *vbox, *hbox, *label;
    GtkWidget *dialog, *image, *button;
    GtkWidget *ebox, *abox;
    gchar *buf, *sysinfo = NULL;

    dialog = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(dialog),_("About gretl")); 
    gtk_window_set_transient_for(GTK_WINDOW(dialog), 
				 GTK_WINDOW(mdata->main));

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    abox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);
    gtk_box_set_spacing(GTK_BOX(vbox), 5);
    gtk_container_set_border_width(GTK_CONTAINER(abox), 5);
    gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_MOUSE);

    /* arrange for a little horizontal padding */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(hbox), 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);
    
    image = open_logo();
    if (image != NULL) {
	gtk_box_pack_start(GTK_BOX(vbox), image, FALSE, FALSE, 20);
	gtk_widget_show(image);
    }

#ifdef PKGBUILD
# if defined(_WIN64)
    sysinfo = g_strdup_printf("MS Windows (x86_64)");
# elif defined(G_OS_WIN32)
    sysinfo = g_strdup_printf("MS Windows (x86)");
# elif defined(MAC_NATIVE)
    sysinfo = g_strdup_printf("Mac OS X (quartz, x86_64)");
# elif defined(__ppc__)
    sysinfo = g_strdup_printf("Mac OS X (X11, ppc)");
# elif defined(OS_OSX)
    sysinfo = g_strdup_printf("Mac OS X (X11, x86)");
# endif
#elif defined(__GNUC__)
# if defined(linux) && defined(__x86_64__)
    sysinfo = g_strdup_printf("Linux x86_64");
# elif defined(__x86_64__)
    sysinfo = g_strdup_printf("Intel x86_64");
# endif
#endif

    if (sysinfo == NULL) {
	buf = g_markup_printf_escaped("<span weight=\"bold\" size=\"x-large\">"
				      "gretl %s</span>\n"
				      "%s %s\n%s", GRETL_VERSION, 
				      _("build date"), BUILD_DATE, 
				      _(bonmot));
    } else {
	buf = g_markup_printf_escaped("<span weight=\"bold\" size=\"x-large\">"
				      "gretl %s</span>\n"
				      "%s\n%s %s\n%s", GRETL_VERSION,
				      sysinfo, _("build date"), BUILD_DATE, 
				      _(bonmot));
	g_free(sysinfo);
    }

    label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label), buf);
    g_free(buf);
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);

    /* Website link */
    ebox = gtk_event_box_new();
    buf = g_markup_printf_escaped("<span color=\"blue\"><u>%s</u></span>",
				  website);
    label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label), buf);
    g_free(buf);
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_container_add(GTK_CONTAINER(ebox), label);
    gtk_box_pack_start(GTK_BOX(vbox), ebox, FALSE, FALSE, 0);
    g_signal_connect(ebox, "button-press-event",
		     G_CALLBACK(show_website), NULL);
    g_signal_connect(ebox, "enter-notify-event",
		     G_CALLBACK(show_link_cursor), NULL);

    /* Copyright label */
    buf = g_markup_printf_escaped("<span size=\"small\">%s</span>",
				  copyright);
    label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label), buf);
    g_free(buf);
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 5);

    /* Translator credits */
    if (strcmp(_("translator_credits"), "translator_credits")) {
	buf = g_markup_printf_escaped("<span size=\"small\">%s</span>",
				      _("translator_credits"));
	label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(label), buf);
	g_free(buf);
	gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
	gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
    }

    /* NEWS button */
    button = gtk_button_new_with_label(_("News"));
    gtk_box_pack_start(GTK_BOX(abox), button, FALSE, FALSE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(relnotes_callback), 
		     NULL);

    /* GPL button */
    button = gtk_button_new_with_label(_("License"));
    gtk_box_pack_start(GTK_BOX(abox), button, FALSE, FALSE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(license_callback), 
		     NULL);

    /* OK button */
    button = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, FALSE, FALSE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(button);

    gtk_widget_show_all(dialog);
}
