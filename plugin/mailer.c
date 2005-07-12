/*
 *  This driver file Copyright (c) 2005 by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* Mailer plugin for gretl, based on mpack, by John G. Myers.
   Please see the files in ./mpack for the Carnegie Mellon
   copyright notice. */

#include "libgretl.h"

#include <gtk/gtk.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>

#include "mpack/mpack.h"

extern int errno;

#if GTK_MAJOR_VERSION < 2

enum {
    GTK_STOCK_OK,
    GTK_STOCK_CANCEL
};

# define G_OBJECT(o)                    GTK_OBJECT(o)
# define g_object_set_data(o,s,d)       gtk_object_set_data(o,s,d)
# define g_object_get_data(o,s)         gtk_object_get_data(o,s)
# define G_CALLBACK(f)                  GTK_SIGNAL_FUNC(f)
# define g_signal_connect(o,s,f,p)      gtk_signal_connect(o,s,f,p)

GtkWidget *standard_button (int code)
{
    const char *button_strings[] = {
	N_("OK"),
	N_("Cancel")
    };

    return gtk_button_new_with_label(_(button_strings[code]));
}

static gint entry_activate (GtkWidget *w, GdkEventKey *key, gpointer p)
{
    GtkWidget *top = gtk_widget_get_toplevel(w);

    gtk_window_activate_default(GTK_WINDOW(top));
    return FALSE;
}

void gtk_entry_set_activates_default (GtkEntry *entry, gboolean setting)
{
    gtk_signal_connect(GTK_OBJECT(entry), "activate", 
		       GTK_SIGNAL_FUNC(entry_activate), NULL);
}

#else

# define standard_button(s) gtk_button_new_from_stock(s)

#endif /* alternate gtk versions */

struct mail_info {
    GtkWidget *dlg;
    GtkWidget *recip_entry;
    GtkWidget *subj_entry;
    GtkWidget *note_entry;
    GtkWidget *ok;
    GtkWidget *cancel;
    char **recip;
    char **subj;
    char **note;
    int *errp;
};

static void finalize_mail_settings (GtkWidget *w, struct mail_info *minfo)
{
    if (w == minfo->cancel) {
	*minfo->errp = -1;
    } else {
	const gchar *txt;
	int err = 0;

	txt = gtk_entry_get_text(GTK_ENTRY(minfo->recip_entry));
	if (txt != NULL && *txt != '\0') {
	    *minfo->recip = g_strdup(txt);
	} else {
	    err = 1;
	}

	if (!err) {
	    txt = gtk_entry_get_text(GTK_ENTRY(minfo->subj_entry));
	    if (txt != NULL && *txt != '\0') {
		*minfo->subj = g_strdup(txt);
	    }

	    txt = gtk_entry_get_text(GTK_ENTRY(minfo->note_entry));
	    if (txt != NULL && *txt != '\0') {
		*minfo->note = g_strdup(txt);
	    }
	}

	*minfo->errp = err;
    }

    gtk_widget_destroy(minfo->dlg);
}

static void set_dialog_border_widths (GtkWidget *dlg)
{
    int w1 = 10, w2 = 5;

#if GTK_MAJOR_VERSION < 2
    gtk_container_border_width(GTK_CONTAINER 
			       (GTK_DIALOG(dlg)->vbox), w1);
    gtk_container_border_width(GTK_CONTAINER 
			       (GTK_DIALOG(dlg)->action_area), w2);
#else
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dlg)->vbox), w1);
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dlg)->action_area), w2);
#endif
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dlg)->vbox), w2);
}

static int 
mail_to_dialog (const char *fname, char **recip, char **subj, char **note)
{
    const gchar *lbls[] = {
	N_("To:"),
	N_("Subject:"),
	N_("Note:")
    };
    GtkWidget *tbl, *lbl, *hbox;
    const char *short_fname, *p;
    gchar *attach_str;
    struct mail_info minfo;
    int i, err = 0;

    minfo.dlg = gtk_dialog_new();
    minfo.recip = recip;
    minfo.subj = subj;
    minfo.note = note;
    minfo.errp = &err;

    g_signal_connect(G_OBJECT(minfo.dlg), "destroy", 
		     G_CALLBACK(gtk_main_quit), NULL);

    gtk_window_set_title(GTK_WINDOW(minfo.dlg), _("gretl: send mail"));
    set_dialog_border_widths(minfo.dlg);
    gtk_window_set_position(GTK_WINDOW(minfo.dlg), GTK_WIN_POS_MOUSE);
    gtk_widget_set_usize(minfo.dlg, 390, 185);

    tbl = gtk_table_new(3, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(minfo.dlg)->vbox), tbl);
   
    for (i=0; i<3; i++) {
	GtkWidget *w;

	lbl = gtk_label_new(_(lbls[i]));
	gtk_table_attach_defaults(GTK_TABLE (tbl), lbl, 0, 1, i, i+1);

	w = gtk_entry_new();
	gtk_entry_set_activates_default(GTK_ENTRY(w), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), w, 1, 2, i, i+1);

	if (i == 0) {
	    minfo.recip_entry = w;
	} else if (i == 1) {
	    minfo.subj_entry = w;
	} else {
	    minfo.note_entry = w;
	}
    }

    short_fname = fname;
    if ((p = strrchr(fname, '/')) != NULL) {
	short_fname = p + 1;
    }

    hbox = gtk_hbox_new(FALSE, 5);
    attach_str = g_strdup_printf(_("sending %s as attachment"), short_fname);
    lbl = gtk_label_new(attach_str);
    g_free(attach_str);
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(minfo.dlg)->vbox), hbox, FALSE, FALSE, 5);

    /* Create the "OK" button */
    minfo.ok = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(minfo.ok, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(minfo.dlg)->action_area), 
		       minfo.ok, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(minfo.ok), "clicked", 
		     G_CALLBACK(finalize_mail_settings), &minfo);
    gtk_widget_grab_default(minfo.ok);

    /* And a Cancel button */
    minfo.cancel = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(minfo.cancel, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(minfo.dlg)->action_area), 
		       minfo.cancel, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(minfo.cancel), "clicked", 
		     G_CALLBACK(finalize_mail_settings), &minfo);

    gtk_widget_show_all(minfo.dlg);
    gtk_window_set_modal(GTK_WINDOW(minfo.dlg), TRUE);
    gtk_main();

    return err;
}

static void sendmail (FILE *infile, char *recipient)
{
    char *args[4];
    int status;
    int pid;

    args[0] = "sendmail";
    args[1] = "-oi";
    args[2] = recipient;
    args[3] = NULL;

    do {
	pid = fork();
    } while (pid == -1 && errno == EAGAIN);
    
    if (pid == -1) {
	perror("fork");
	return;
    }

    if (pid != 0) {
	while (pid != wait(&status));
	return;
    }

    dup2(fileno(infile), 0);
    fclose(infile);
    execv("/usr/sbin/sendmail", args);
    execv("/usr/lib/sendmail", args);
    perror("execv");
    _exit(1);
}

static int pack_and_mail (const char *subject, const char *note,
			  const char *fname, const char *ctype,
			  char *recipient, char *tmpfname)
{
    FILE *fp;
    int err = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	perror(fname);
	err = 1;
    }

    if (!err) {
	err = encode(fp, fname, note, subject, recipient,
		     ctype, tmpfname);
    }

    if (!err) {
	fp = fopen(tmpfname, "r");
	if (fp == NULL) {
	    perror(tmpfname);
	    err = 1;
	} else {
	    sendmail(fp, recipient);
	    fclose(fp);
	}
    }

    remove(tmpfname);

    return err;
}

int email_file (const char *fname, const char *userdir, char *errmsg)
{
    char temp[FILENAME_MAX];
    char *subject = NULL;
    char *note = NULL;
    char *recipient = NULL;
    char *ctype = "application/x-gretldata";
    int err = 0;

    *errmsg = 0;

    sprintf(temp, "%smpack.XXXXXX", userdir);
    if (mktemp(temp) == NULL) {
	err = 1;
    }

    if (!err) {
	err = mail_to_dialog(fname, &recipient, &subject, &note);
    }

    if (!err) {
	err = pack_and_mail(subject, note, fname, ctype, recipient,
			    temp);
    }

    return err;
}

