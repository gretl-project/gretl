/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */


/* The "about" dialogs for the several gretl GUI variants. */

#include "gretl.h"
#include "version.h"

#ifdef G_OS_WIN32 
# include "build.h"
#endif

/* some material that is common between gtk 1.2 and gtk 2 variants */

const gchar *copyright = "Copyright (C) 2000-2005 Allin Cottrell";
const gchar *website = "http://gretl.sourceforge.net/";

#ifdef USE_GNOME
const gchar *
gretl_gnome_blurb = N_("An econometrics program for the gnome desktop "
		       "issued under the GNU General Public License.  "
		       "http://gretl.sourceforge.net/");
#endif

/* end of common stuff */

#ifndef OLD_GTK  /* we'll start with the current versions */

# ifdef USE_GNOME

void about_dialog (gpointer data)
{
    static GtkWidget *about = NULL;
    gchar *pixfile;
    GdkPixbuf* pbuf = NULL;
	
    gchar *authors[] = {
	"Allin Cottrell <cottrell@wfu.edu>",
	NULL
    };
    gchar *documenters[] = {
	"Allin Cottrell <cottrell@wfu.edu>",
	NULL
    };
    gchar *translator_credits = _("translator_credits");

    if (about != NULL) {
	gdk_window_show (about->window);
	gdk_window_raise (about->window);
	return;
    }

    pixfile = gnome_program_locate_file(NULL,
					GNOME_FILE_DOMAIN_PIXMAP,
					"gretl-logo.xpm",
					TRUE,
					NULL);

    if (pixfile != NULL) {
	pbuf = gdk_pixbuf_new_from_file(pixfile, NULL);
    } else {
	fprintf(stderr, "Couldn't find gretl-logo.xpm\n");
    }

    about = gnome_about_new ("gretl", GRETL_VERSION,
			     copyright,
			     _(gretl_gnome_blurb),
			     (const char **)authors,
			     (const char **)documenters,
			     strcmp (translator_credits, "translator_credits") != 0 ?
			     (const char *)translator_credits : NULL,
			     pbuf);

    gtk_window_set_transient_for (GTK_WINDOW (about),
				  GTK_WINDOW (mdata->w));

    gtk_window_set_destroy_with_parent (GTK_WINDOW (about), TRUE);

    if (pbuf != NULL)
	g_object_unref(pbuf);
	
    g_signal_connect (G_OBJECT (about), "destroy",
		      G_CALLBACK (gtk_widget_destroyed), &about);
	
    gtk_widget_show (about);
}

# else /* plain GTK 2 version of About dialog follows */

static GtkWidget *open_logo (const char *pngname)
{
    char fullname[MAXLEN];
    GdkPixbuf *pbuf;
    GError *error = NULL;
    GtkWidget *image;

    build_path(paths.gretldir, pngname, fullname, NULL);

    pbuf = gdk_pixbuf_new_from_file (fullname, &error);

    if (pbuf == NULL) {
	errbox(error->message);
	g_error_free(error);
	return NULL;
    } else {
	image = gtk_image_new_from_pixbuf (pbuf);
	return image;
    }
}

static void about_table_setup (GtkWidget *vbox, GtkWidget *view)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new (NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vbox), 
                       sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (sw),
                                    GTK_POLICY_AUTOMATIC,
                                    GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW (sw),
                                         GTK_SHADOW_IN);
    gtk_container_add (GTK_CONTAINER(sw), view); 
    gtk_widget_show(view);
    gtk_widget_show(sw);
}

void about_dialog (gpointer data) 
{
    GtkWidget *notebook, *box, *label, *tempwid;
    GtkWidget *view, *dialog;
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    char *tempstr, *no_gpl, buf[MAXSTR];
    const gchar *tr_credit = "";
    FILE *fd;

    no_gpl = 
	g_strdup_printf (_("Cannot find the license agreement file COPYING. "
			   "Please make sure it's in %s"), 
			 paths.gretldir);
    dialog = gtk_dialog_new ();
    gtk_window_set_title(GTK_WINDOW(dialog),_("About gretl")); 
    gtk_container_set_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
      
    notebook = gtk_notebook_new ();
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			notebook, TRUE, TRUE, 0);

    /* construct the first page */
    box = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    if ((tempwid = open_logo("gretl-logo.xpm"))) {
	gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 30);
	gtk_widget_show (tempwid);
    }

#  ifdef ENABLE_NLS
    if (strcmp(_("translator_credits"), "translator_credits")) {
	tr_credit = _("translator_credits");
    }
#  endif    
    
    tempstr = g_strdup_printf ("gretl, version %s\n"
#  ifdef G_OS_WIN32
			       BUILD_DATE
#  endif
			       "%s\n%s\n%s",
			       GRETL_VERSION, copyright, 
			       website, tr_credit);
    tempwid = gtk_label_new (tempstr);
    g_free (tempstr);

    gtk_label_set_justify(GTK_LABEL(tempwid), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 0);
    gtk_widget_show (tempwid);

    gtk_widget_show(box);

    label = gtk_label_new (_("About"));
    gtk_widget_show (label);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    /* now the second page */
    box = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (box), 10);

    view = gtk_text_view_new ();
    gtk_text_view_set_editable (GTK_TEXT_VIEW (view), FALSE);
    gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (view), GTK_WRAP_NONE);
    gtk_widget_modify_font(GTK_WIDGET(view), fixed_font);

    about_table_setup(box, view);

    gtk_widget_show (box);

    label = gtk_label_new (_("License Agreement"));
    gtk_widget_show (label);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    tempwid = gtk_button_new_from_stock(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			tempwid, FALSE, FALSE, 0);
    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      G_CALLBACK (delete_widget), 
		      dialog);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_iter_at_offset (tbuf, &iter, 0);

    tempstr = g_strdup_printf("%s/COPYING", paths.gretldir);
    if ((fd = gretl_fopen(tempstr, "r")) == NULL) {
	gtk_text_buffer_insert(tbuf, &iter, no_gpl, -1);
	gtk_widget_show(dialog);
	g_free(tempstr);
	return;
    }
    g_free(tempstr);
   
    memset(buf, 0, sizeof (buf));
    while(fread (buf, 1, sizeof (buf) - 1, fd)) {
	gtk_text_buffer_insert(tbuf, &iter, buf, strlen (buf));
	memset(buf, 0, sizeof (buf));
    }
    fclose (fd);

    gtk_widget_show(notebook);
    gtk_widget_set_size_request(dialog, 520, 420);
    gtk_widget_show(dialog);

    g_free(no_gpl);
}
         
# endif /* end of gnome/plain gtk variants for gtk 2 */

#else /* now on to the old gtk variants */

# ifdef USE_GNOME

void about_dialog (gpointer data) 
{
    GtkWidget* dlg;
    char const *authors[] = {
	"Allin Cottrell",
	NULL
    };
    gchar *comment = NULL;

#  ifdef ENABLE_NLS
    if (strcmp(_("translator_credits"), "translator_credits")) {
	comment = g_strconcat(_(gretl_gnome_blurb), " ", _("translator_credits"),
			      NULL);
    }
#  endif 

    dlg = gnome_about_new("gretl", GRETL_VERSION,
			  copyright, 
			  authors, 
			  (comment != NULL)? comment : _(gretl_gnome_blurb),
			  gnome_pixmap_file("gretl-logo.xpm") 
			  );

    if (comment != NULL) g_free(comment);

    gnome_dialog_set_parent(GNOME_DIALOG(dlg), GTK_WINDOW(mdata->w));

    gtk_widget_show(dlg);
}

# else /* plain gtk 1.2 version of About dialog follows */

static int open_xpm (char *filename, GtkWidget *parent, GdkPixmap **pixmap, 
		     GdkBitmap **mask) 
{
    char exfile[MAXLEN];
    GtkStyle *style;

    if (*filename == '\0') return 1;
    strcpy(exfile, paths.gretldir);
    if (exfile[strlen(exfile) - 2] != SLASH)
	strcat(exfile, SLASHSTR);
    strcat(exfile, filename);

    style = gtk_widget_get_style (parent);
    *pixmap = gdk_pixmap_create_from_xpm (parent->window, 
					  mask, 
					  &style->bg[GTK_STATE_NORMAL], 
					  exfile);

    return (*pixmap != NULL);
}

void about_dialog (gpointer data) 
{
    GtkWidget *tempwid, *notebook, *box, *label, *view, *vscroll;
    GdkPixmap *logo_pixmap;
    GdkBitmap *logo_mask;
    char *tempstr, *no_gpl, buf[MAXSTR];
    const gchar *tr_credit = "";
    GtkWidget *dialog;
    FILE *fd;

    no_gpl = 
	g_strdup_printf (_("Cannot find the license agreement file COPYING. "
			 "Please make sure it's in %s"), 
			 paths.gretldir);
    dialog = gtk_dialog_new ();
    gtk_window_set_title (GTK_WINDOW (dialog), _("About gretl"));
    gtk_widget_set_usize(dialog, 522, 470);

    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
    gtk_box_set_homogeneous (GTK_BOX (GTK_DIALOG (dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
    gtk_signal_connect_object (GTK_OBJECT (dialog), 
			       "delete_event", GTK_SIGNAL_FUNC 
			       (gtk_widget_destroy), GTK_OBJECT (dialog));
    gtk_widget_realize (dialog);
      
    notebook = gtk_notebook_new ();
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			notebook, TRUE, TRUE, 0);
    gtk_widget_show (notebook);
   
    box = gtk_vbox_new (TRUE, 5);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);
   
    if (open_xpm ("gretl-logo.xpm", mdata->w, &logo_pixmap, &logo_mask)) {
	tempwid = gtk_pixmap_new (logo_pixmap, logo_mask);
	gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 0);
	gtk_widget_show (tempwid);
    }

#  ifdef ENABLE_NLS
    if (strcmp(_("translator_credits"), "translator_credits")) {
	tr_credit = _("translator_credits");
    }
#  endif  

    tempstr = g_strdup_printf ("gretl, version %s\n"
			       "%s\n%s\n%s",
			       GRETL_VERSION, copyright, 
			       website, tr_credit);
    tempwid = gtk_label_new (tempstr);
    g_free (tempstr);

    gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 0);
    gtk_widget_show (tempwid);
   
    label = gtk_label_new (_("About"));
    gtk_widget_show (label);
   
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    box = gtk_vbox_new (FALSE, 5);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    tempwid = gtk_table_new (1, 2, FALSE);
    gtk_box_pack_start (GTK_BOX (box), tempwid, TRUE, TRUE, 0);
    gtk_widget_show (tempwid);

    view = gtk_text_new (NULL, NULL);
    gtk_text_set_editable (GTK_TEXT (view), FALSE);
    gtk_text_set_word_wrap (GTK_TEXT (view), TRUE);
    gtk_table_attach (GTK_TABLE (tempwid), view, 0, 1, 0, 1,
		      GTK_FILL | GTK_EXPAND, GTK_FILL | 
		      GTK_EXPAND | GTK_SHRINK, 0, 0);
    gtk_widget_show (view);

    vscroll = gtk_vscrollbar_new (GTK_TEXT (view)->vadj);
    gtk_table_attach (GTK_TABLE (tempwid), vscroll, 1, 2, 0, 1,
		      GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);
    gtk_widget_show (vscroll);

    label = gtk_label_new (_("License Agreement"));
    gtk_widget_show(label);
   
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    tempwid = gtk_button_new_with_label (_("  Close  "));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			tempwid, FALSE, FALSE, 0);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    tempstr = g_strdup_printf("%s/COPYING", paths.gretldir);
    if ((fd = gretl_fopen(tempstr, "r")) == NULL) {
	gtk_text_insert(GTK_TEXT(view), NULL, NULL, NULL, 
			no_gpl, strlen(no_gpl));
	gtk_widget_show(dialog);
	g_free(tempstr);
	return;
    }
    g_free(tempstr);
   
    memset(buf, 0, sizeof buf);
    while (fread(buf, 1, sizeof buf - 1, fd)) {
	gtk_text_insert(GTK_TEXT (view), 
			fixed_font, NULL, NULL, buf, strlen(buf));
	memset (buf, 0, sizeof buf);
    }
    fclose(fd);

    gtk_widget_show(dialog);

    g_free(no_gpl);
} 
        
# endif /* not GNOME */

#endif /* end of gtk version alternation */
