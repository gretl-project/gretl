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
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "mpack/mpack.h"

enum {
    MAIL_OK,
    MAIL_NO_RECIPIENT,
    MAIL_CANCEL
};

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
    GtkWidget *recip_combo;
    GtkWidget *reply_entry;
    GtkWidget *subj_entry;
    GtkWidget *note_entry;
    GtkWidget *ok;
    GtkWidget *cancel;
    char **recip;
    char **reply_to;
    char **subj;
    char **note;
    int *errp;
    char *addr_file;
    GList *addrs;
    char *orig_reply_to;
};

static void destroy_mail_info (GtkWidget *w, struct mail_info *minfo)
{
    g_free(minfo->addr_file);
    g_free(minfo->orig_reply_to);
    /* and the GList?? */

    gtk_main_quit();
}

static void save_gretl_addresses (struct mail_info *minfo)
{
    FILE *fp;

    fp = gretl_fopen(minfo->addr_file, "w");

    if (fp != NULL) {
	GList *list = minfo->addrs;
	int i, maxaddrs = 10;

	if (minfo->reply_to != NULL && *minfo->reply_to != NULL &&
	    **minfo->reply_to != '\0') {
	    fprintf(fp, "Reply-To: %s\n", *minfo->reply_to);
	}
	for (i=0; i<maxaddrs && list != NULL; i++) {
	    fprintf(fp, "%s\n", (char *) list->data);
	    list = list->next;
	}
	fclose(fp);
    } 
}

static void finalize_mail_settings (GtkWidget *w, struct mail_info *minfo)
{
    int save = 0;

    if (w == minfo->cancel) {
	*minfo->errp = MAIL_CANCEL;
    } else {
	GList *list = minfo->addrs;
	const gchar *txt;
	int err = MAIL_OK;

	txt = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(minfo->recip_combo)->entry));
	if (txt != NULL && *txt != '\0') {
	    *minfo->recip = g_strdup(txt);
	    save = 1;
	    while (list) {
		if (!strcmp(txt, (char *) list->data)) {
		    save = 0;
		    break;
		}
		list = g_list_next(list);
	    }
	    if (save) {
		minfo->addrs = g_list_prepend(minfo->addrs, g_strdup(txt));
	    }
	} else {
	    err = MAIL_NO_RECIPIENT;
	}

	if (!err) {
	    txt = gtk_entry_get_text(GTK_ENTRY(minfo->reply_entry));
	    if (txt != NULL && *txt != '\0') {
		*minfo->reply_to = g_strdup(txt);
		if (minfo->orig_reply_to == NULL ||
		    strcmp(txt, minfo->orig_reply_to)) {
		    save = 1;
		}
	    }

	    txt = gtk_entry_get_text(GTK_ENTRY(minfo->subj_entry));
	    if (txt != NULL && *txt != '\0') {
		*minfo->subj = g_strdup(txt);
	    }

	    txt = gtk_entry_get_text(GTK_ENTRY(minfo->note_entry));
	    if (txt != NULL && *txt != '\0') {
		*minfo->note = g_strdup_printf("%s\n", txt);
	    }
	}

	*minfo->errp = err;
    }

    if (save) {
	save_gretl_addresses(minfo);
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

static void get_gretl_address_list (struct mail_info *minfo)
{
    GList *addrs = NULL;
    gchar *reply = NULL;
    FILE *fp;

    minfo->addr_file = g_strdup_printf("%sgretl.addresses", gretl_user_dir());

    fp = gretl_fopen(minfo->addr_file, "r");
    if (fp != NULL) {
	char line[128];

	while (fgets(line, sizeof line, fp)) {
	    if (string_is_blank(line)) {
		continue;
	    }
	    line[strlen(line) - 1] = '\0';
	    if (!strncmp(line, "Reply-To:", 9)) {
		reply = g_strdup(line + 10);
	    } else {
		addrs = g_list_append(addrs, g_strdup(line));
	    }
	}

	fclose(fp);
    } 

    minfo->addrs = addrs;
    minfo->orig_reply_to = reply;
}

static void cancel_mail (struct mail_info *minfo)
{
    fprintf(stderr, "delete-event: canceling\n");
    *minfo->errp = MAIL_CANCEL;
}

static int mail_to_dialog (const char *fname, char **recip, char **reply_to, 
			   char **subj, char **note)
{
    const gchar *lbls[] = {
	N_("To:"),
	N_("Reply-To:"),
	N_("Subject:"),
	N_("Note:")
    };
    GtkWidget *tbl, *lbl, *hbox;
    const char *short_fname, *p;
    gchar *attach_str;
    struct mail_info minfo;
    int i, err = 0;

    get_gretl_address_list(&minfo);

    minfo.dlg = gtk_dialog_new();

    minfo.recip = recip;
    minfo.reply_to = reply_to;
    minfo.subj = subj;
    minfo.note = note;
    minfo.errp = &err;

    g_signal_connect(G_OBJECT(minfo.dlg), "delete-event", 
		     G_CALLBACK(cancel_mail), &minfo);

    g_signal_connect(G_OBJECT(minfo.dlg), "destroy", 
		     G_CALLBACK(destroy_mail_info), &minfo);

    gtk_window_set_title(GTK_WINDOW(minfo.dlg), _("gretl: send mail"));
    set_dialog_border_widths(minfo.dlg);
    gtk_window_set_position(GTK_WINDOW(minfo.dlg), GTK_WIN_POS_MOUSE);
    gtk_widget_set_usize(minfo.dlg, 420, 224);

    tbl = gtk_table_new(4, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(minfo.dlg)->vbox), tbl);

    short_fname = fname;
    if ((p = strrchr(fname, '/')) != NULL) {
	short_fname = p + 1;
    }    

    for (i=0; i<4; i++) {
	GtkWidget *w;

	lbl = gtk_label_new(_(lbls[i]));
	gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
	gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, i, i+1, GTK_FILL, GTK_FILL, 0, 0);

	if (i == 0) {
	    w = gtk_combo_new();
	    if (minfo.addrs != NULL) {
		gtk_combo_set_popdown_strings(GTK_COMBO(w), minfo.addrs);
	    } 
	} else {
	    w = gtk_entry_new();
	}

	if (i == 1) {
	    if (minfo.orig_reply_to != NULL) {
		gtk_entry_set_text(GTK_ENTRY(w), minfo.orig_reply_to);
	    }
	} else if (i == 2) {
	    gtk_entry_set_text(GTK_ENTRY(w), "dataset");
	} else if (i == 3) {
	    gchar *note = g_strdup_printf("Please find the gretl data file %s attached.",
					  short_fname);

	    gtk_entry_set_text(GTK_ENTRY(w), note);
	    g_free(note);
	}

	if (i == 0) {
	    gtk_entry_set_activates_default(GTK_ENTRY(GTK_COMBO(w)->entry), TRUE);
	} else {
	    gtk_entry_set_activates_default(GTK_ENTRY(w), TRUE);
	}

	gtk_table_attach_defaults(GTK_TABLE(tbl), w, 1, 2, i, i+1);

	if (i == 0) {
	    minfo.recip_combo = w;
	} else if (i == 1) {
	    minfo.reply_entry = w;
	} else if (i == 2) {
	    minfo.subj_entry = w;
	} else {
	    minfo.note_entry = w;
	}
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

#if GTK_MAJOR_VERSION < 2

static void mail_infobox (const char *msg) 
{
    GtkWidget *w, *label, *button, *vbox, *hbox;

    w = gtk_window_new(GTK_WINDOW_DIALOG);

    gtk_container_border_width(GTK_CONTAINER(w), 5);
    gtk_window_position (GTK_WINDOW(w), GTK_WIN_POS_MOUSE);
    gtk_window_set_title (GTK_WINDOW (w), _("gretl info")); 

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(w), vbox);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);

    /* text of message */
    label = gtk_label_new(msg);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    /* button */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    
    button = gtk_button_new_with_label(_("OK"));

    gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 5);

    gtk_signal_connect_object(GTK_OBJECT(button), "clicked",
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      (gpointer) w);

    gtk_widget_show_all(w);
}

#else

static void mail_infobox (const char *msg)
{
    GtkWidget *dialog;

    dialog = gtk_message_dialog_new (NULL,
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     GTK_MESSAGE_INFO,
				     GTK_BUTTONS_CLOSE,
				     msg);
    gtk_dialog_run(GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
}

#endif


static void
real_send_mail (FILE *infile, char *recipient, char *sendmail)
{
    char *argv[4];
    int status;
    int pid;

    extern int errno;

    argv[0] = "sendmail";
    argv[1] = "-oi";
    argv[2] = recipient;
    argv[3] = NULL;

    do {
	pid = fork();
    } while (pid == -1 && errno == EAGAIN);
    
    if (pid == -1) {
	perror("fork");
	return;
    }

    if (pid != 0) {
	/* parent */
	while (pid != wait(&status));
	if (WIFEXITED(status)) {
	    int err = WEXITSTATUS(status);
	    gchar *msg = g_strdup_printf("sendmail exited with status %d", err);

	    mail_infobox(msg);
	    g_free(msg);
	} else {
	    mail_infobox("sendmail exited abnormally");
	}
	return;
    }

    /* child */
#if 0
    fprintf(stderr, "not doing sendmail\n");
    fclose(infile);
#else
    dup2(fileno(infile), 0);
    fclose(infile);
    execv(sendmail, argv);
    perror("execv");
#endif
    _exit(1);
}

static int pack_and_mail (const char *subject, const char *note,
			  const char *fname, const char *ctype,
			  char *recipient, char *reply_to,
			  char *tmpfname, char *sendmail)
{
    FILE *fp;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	perror(fname);
	err = 1;
    }

    if (!err) {
	err = encode(fp, fname, note, subject, recipient,
		     reply_to, ctype, tmpfname);
    }

    if (!err) {
	fp = gretl_fopen(tmpfname, "r");
	if (fp == NULL) {
	    perror(tmpfname);
	    err = 1;
	} else {
	    real_send_mail(fp, recipient, sendmail);
	    fclose(fp);
	}
    }

    remove(tmpfname);

    return err;
}

static char *find_sendmail (void)
{
    const char *candidate[] = {
	"/usr/sbin/sendmail",
	"/usr/lib/sendmail",
	"/usr/bin/sendmail",
	NULL
    };
    struct stat buf;
    char *prog = NULL;
    int i;

    for (i=0; candidate[i] != NULL; i++) {
	if (stat(candidate[i], &buf) == 0 && (buf.st_mode & S_IXUSR)) {
	    prog = g_strdup(candidate[i]);
	    break;
	}
    }

    return prog;
}

int email_file (const char *fname, const char *userdir, char *errmsg)
{
    char temp[FILENAME_MAX];
    char *subject = NULL;
    char *note = NULL;
    char *recipient = NULL;
    char *reply_to = NULL;
    char *sendmail = NULL;
    char *ctype = "application/x-gretldata";
    int mval, err = 0;

    *errmsg = 0;

    sendmail = find_sendmail();
    if (sendmail == NULL) {
	strcpy(errmsg, "Couldn't find sendmail executable");
	return 1;
    }

    sprintf(temp, "%smpack.XXXXXX", userdir);
    if (mktemp(temp) == NULL) {
	err = 1;
    }

    if (!err) {
	mval = mail_to_dialog(fname, &recipient, &reply_to, &subject, &note);
	if (mval == MAIL_NO_RECIPIENT) {
	    strcpy(errmsg, "No address was given");
	    err = 1;
	} else if (mval == MAIL_OK) {
	    err = pack_and_mail(subject, note, fname, ctype, recipient,
				reply_to, temp, sendmail);
	}
    }

    g_free(reply_to);
    g_free(sendmail);

    return err;
}

