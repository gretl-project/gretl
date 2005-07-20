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
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>

extern int h_errno;

#include "mpack/mpack.h"

enum {
    MAIL_OK,
    MAIL_NO_RECIPIENT,
    MAIL_NO_SERVER,
    MAIL_NO_SENDER,
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
# define gtk_widget_set_size_request(g,w,h) gtk_widget_set_usize(g,w,h)
# define gtk_notebook_set_current_page(n,p) gtk_notebook_set_page(n,p)

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

struct msg_info {
    char *recip;
    char *sender;
    char *subj;
    char *note;  
};  

struct mail_info {
    char *sender;
    char *sig;
    int want_sig;
    char *server;
    unsigned short port;
};

struct mail_dialog {
    GtkWidget *dlg;
    GtkWidget *recip_combo;
    GtkWidget *reply_entry;
    GtkWidget *subj_entry;
    GtkWidget *note_entry;
    GtkWidget *server_entry;
    GtkWidget *port_entry;
    GtkWidget *ok;
    GtkWidget *cancel;
    struct mail_info *minfo;
    struct msg_info *msg;
    char *addr_file;
    GList *addrs;
    int *errp;
};

static void msg_init (struct msg_info *msg)
{
    msg->recip = NULL;
    msg->sender = NULL;
    msg->subj = NULL;
    msg->note = NULL;
}

static void free_msg (struct msg_info *msg)
{
    free(msg->recip);
    free(msg->sender);
    free(msg->subj);
    free(msg->note);
}

static void mail_info_init (struct mail_info *minfo)
{
    minfo->sender = NULL;
    minfo->sig = NULL;
    minfo->want_sig = 1;
    minfo->server = NULL;
    minfo->port = 25;
}

static void free_mail_info (struct mail_info *minfo)
{
    free(minfo->sender);
    free(minfo->sig);
    free(minfo->server);
}

static void save_email_info (struct mail_dialog *md)
{
    struct mail_info *minfo = md->minfo;
    FILE *fp;

    fp = gretl_fopen(md->addr_file, "w");

    if (fp != NULL) {
	GList *list = md->addrs;
	int i, maxaddrs = 10;

	if (minfo->sender != NULL && *minfo->sender != '\0') {
	    fprintf(fp, "Reply-To: %s\n", minfo->sender);
	}

	if (minfo->server != NULL && *minfo->server != '\0') {
	    fprintf(fp, "SMTP server: %s\n", minfo->server);
	}

	if (minfo->port != 25) {
	    fprintf(fp, "SMTP port: %d\n", minfo->port);
	}

	for (i=0; i<maxaddrs && list != NULL; i++) {
	    fprintf(fp, "%s\n", (char *) list->data);
	    list = list->next;
	}

	fclose(fp);
    } 
}

static char *add_to_string (char *str, const char *add)
{
    char *tmp = NULL;

    if (str == NULL) {
	tmp = malloc(strlen(add) + 1);
	if (tmp != NULL) {
	    *tmp = '\0';
	}
    } else {
	tmp = realloc(str, strlen(str) + strlen(add) + 1);
    }

    if (tmp != NULL) {
	str = tmp;
	strcat(str, add);
    }

    return str;
}

static char *get_signature (void)
{
    char *home = getenv("HOME");
    FILE *fp;
    char *sig = NULL;

    if (home != NULL) {
	gchar *sigfile = g_strdup_printf("%s/.signature", home);
	
	fp = gretl_fopen(sigfile, "r");
	if (fp != NULL) {
	    char line[128];

	    while (fgets(line, sizeof line, fp)) {
		sig = add_to_string(sig, line);
	    }
	    fclose(fp);
	}
	g_free(sigfile);
    }

    return sig;
}

static void finalize_mail_settings (GtkWidget *w, struct mail_dialog *md)
{
    GList *list = md->addrs;
    struct mail_info *minfo = md->minfo;
    struct msg_info *msg = md->msg;
    const gchar *txt;
    int err = MAIL_OK;
    int save = 0;

    if (w == md->cancel) {
	*md->errp = MAIL_CANCEL;
	gtk_widget_destroy(md->dlg);
	return;
    }

    /* recipient */
    txt = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(md->recip_combo)->entry));
    if (txt != NULL && *txt != '\0') {
	int i = 0;

	msg->recip = g_strdup(txt);
	fprintf(stderr, "targ = '%s'\n", msg->recip);
	save = 1;
	while (list) {
	    if (!strcmp(txt, (char *) list->data)) {
		if (i == 0) {
		    /* current recipient is top of the list already */
		    save = 0;
		} else {
		    /* current recipient should be moved to top */
		    list = g_list_remove(list, list->data);
		}
		break;
	    }
	    list = g_list_next(list);
	    i++;
	}
	if (save) {
	    md->addrs = g_list_prepend(md->addrs, g_strdup(txt));
	} 
    } else {
	err = MAIL_NO_RECIPIENT;
    }

    if (!err) {
	/* reply-to address */
	txt = gtk_entry_get_text(GTK_ENTRY(md->reply_entry));
	if (txt != NULL && *txt != '\0') {
	    msg->sender = g_strdup(txt);
	    if (minfo->sender == NULL ||
		strcmp(txt, minfo->sender)) {
		save = 1;
	    }
	    fprintf(stderr, "sender = '%s'\n", msg->sender);
	} else {
	    err = MAIL_NO_SENDER;
	}
    }

    if (!err) {
	/* message subject */
	txt = gtk_entry_get_text(GTK_ENTRY(md->subj_entry));
	if (txt != NULL && *txt != '\0') {
	    msg->subj = g_strdup(txt);
	    fprintf(stderr, "subj = '%s'\n", msg->subj);
	}

	/* message text */
	txt = gtk_entry_get_text(GTK_ENTRY(md->note_entry));
	if (txt != NULL && *txt != '\0') {
	    if (minfo->sig != NULL && !minfo->want_sig) {
		free(minfo->sig);
		minfo->sig = NULL;
	    }
	    if (minfo->sig != NULL) {
		msg->note = g_strdup_printf("%s\n--\n%s\n", txt, minfo->sig);
	    } else {
		msg->note = g_strdup_printf("%s\n", txt);
	    }
	}
 
	/* SMTP server */
	txt = gtk_entry_get_text(GTK_ENTRY(md->server_entry));
	if (txt != NULL && *txt != '\0') {
	    minfo->server = g_strdup(txt);
	    save = 1;
	    fprintf(stderr, "server = '%s'\n", minfo->server);
	} else {
	    err = MAIL_NO_SERVER;
	}
    }

    if (!err) {
	/* port number */
	txt = gtk_entry_get_text(GTK_ENTRY(md->port_entry));
	if (txt != NULL && *txt != '\0') {
	    minfo->port = atoi(txt);
	    if (minfo->port != 25) {
		save = 1;
	    }
	}	    
    }

    *md->errp = err;

    if (save) {
	save_email_info(md);
    }

    gtk_widget_destroy(md->dlg);
}

static void border_width (GtkWidget *w, int b)
{
#if GTK_MAJOR_VERSION < 2
    gtk_container_border_width(GTK_CONTAINER(w), b);
#else
    gtk_container_set_border_width(GTK_CONTAINER(w), b); 
#endif
}

static void set_dialog_border_widths (GtkWidget *dlg)
{
    int w1 = 10, w2 = 5;

    border_width(GTK_DIALOG(dlg)->vbox, w1);
    border_width(GTK_DIALOG(dlg)->action_area, w2);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dlg)->vbox), w2);
}

static void get_email_info (struct mail_dialog *md)
{
    struct mail_info *minfo = md->minfo;
    GList *addrs = NULL;
    FILE *fp;

    md->addr_file = g_strdup_printf("%sgretl.addresses", gretl_user_dir());

    fp = gretl_fopen(md->addr_file, "r");
    if (fp != NULL) {
	char line[128];

	while (fgets(line, sizeof line, fp)) {
	    if (string_is_blank(line)) {
		continue;
	    }
	    chopstr(line);
	    if (!strncmp(line, "Reply-To:", 9)) {
		minfo->sender = g_strdup(line + 10);
	    } else if (!strncmp(line, "SMTP server:", 12)) {
		minfo->server = g_strdup(line + 13);
	    } else if (!strncmp(line, "SMTP port:", 10)) {
		minfo->port = atoi(line + 11);
	    } else {
		addrs = g_list_append(addrs, g_strdup(line));
	    }
	}

	fclose(fp);
    } 

    md->addrs = addrs;
}

static void cancel_mail (struct mail_dialog *md)
{
    fprintf(stderr, "delete-event: canceling\n");
    *md->errp = MAIL_CANCEL;
}

static int is_data_file (const char *fname)
{
    int ret = 1;

    if (fname != NULL && strlen(fname) > 4) {
	ret = !strcmp(fname + strlen(fname) - 4, ".gdt");
    }

    return ret;
}

static void sig_callback (GtkWidget *w, struct mail_info *minfo)
{
    minfo->want_sig = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
}

static void destroy_mail_dialog (GtkWidget *w, struct mail_dialog *md)
{
    GList *tmp = md->addrs;

    while (tmp) {
	/* free the address strings */
	g_free(tmp->data);
	tmp = g_list_next(tmp);
    }

    g_free(md->addr_file);

    gtk_main_quit();
}

static int 
mail_to_dialog (const char *fname, struct mail_info *minfo, struct msg_info *msg)
{
    const gchar *lbls[] = {
	N_("To:"),
	N_("Reply-To:"),
	N_("Subject:"),
	N_("Note:")
    };
    GtkWidget *tbl, *lbl, *vbox;
    GtkWidget *nb, *hbox;
    gchar *port_str;
    const char *short_fname, *p;
    struct mail_dialog md;
    int datafile, nrows;
    int i, err = 0;

    md.dlg = gtk_dialog_new();
    md.minfo = minfo;
    md.msg = msg;
    md.errp = &err;

    get_email_info(&md);
    md.minfo->sig = get_signature();
    md.minfo->want_sig = minfo->sig != NULL;

    g_signal_connect(G_OBJECT(md.dlg), "delete-event", 
		     G_CALLBACK(cancel_mail), &md);

    g_signal_connect(G_OBJECT(md.dlg), "destroy", 
		     G_CALLBACK(destroy_mail_dialog), &md);

    gtk_window_set_title(GTK_WINDOW(md.dlg), _("gretl: send mail"));
    set_dialog_border_widths(md.dlg);
    gtk_window_set_position(GTK_WINDOW(md.dlg), GTK_WIN_POS_MOUSE);

    nb = gtk_notebook_new();
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(md.dlg)->vbox), nb);
    hbox = gtk_hbox_new(FALSE, 5);
    border_width(hbox, 5);
    vbox = gtk_vbox_new(FALSE, 5);
    border_width(vbox, 5);
    gtk_container_add(GTK_CONTAINER(hbox), vbox);
    lbl = gtk_label_new(_("Message"));
    gtk_notebook_append_page(GTK_NOTEBOOK(nb), hbox, lbl);    

    nrows = (md.minfo->sig == NULL)? 4 : 5;

    tbl = gtk_table_new(nrows, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_container_add(GTK_CONTAINER(vbox), tbl);

    short_fname = fname;
    if ((p = strrchr(fname, '/')) != NULL) {
	short_fname = p + 1;
    }  

    datafile = is_data_file(short_fname);

    for (i=0; i<4; i++) {
	GtkWidget *w;

	lbl = gtk_label_new(_(lbls[i]));
	gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
	gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, i, i+1, GTK_FILL, GTK_FILL, 0, 0);

	if (i == 0) {
	    w = gtk_combo_new();
	    if (md.addrs != NULL) {
		gtk_combo_set_popdown_strings(GTK_COMBO(w), md.addrs);
	    } 
	} else {
	    w = gtk_entry_new();
	}

	if (i == 1) {
	    if (md.minfo->sender != NULL) {
		gtk_entry_set_text(GTK_ENTRY(w), md.minfo->sender);
	    }
	} else if (i == 2) {
	    gtk_entry_set_text(GTK_ENTRY(w), (datafile)? "dataset" : "script");
	} else if (i == 3) {
	    gchar *note;

	    if (datafile) {
		note = g_strdup_printf("Please find the gretl data file %s attached.",
				       short_fname);
	    } else {
		note = g_strdup_printf("Please find the gretl script %s attached.",
				       short_fname);
	    }		
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
	    md.recip_combo = w;
	} else if (i == 1) {
	    md.reply_entry = w;
	} else if (i == 2) {
	    md.subj_entry = w;
	} else {
	    md.note_entry = w;
	}
    }

    if (md.minfo->sig != NULL) {
	GtkWidget *w;

	w = gtk_check_button_new_with_label("Append signature");
	g_signal_connect(G_OBJECT(w), "toggled", G_CALLBACK(sig_callback), minfo);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), w, 0, 2, 4, 5);
	i++;
    }

    hbox = gtk_hbox_new(FALSE, 5);
    border_width(hbox, 5);
    vbox = gtk_vbox_new(FALSE, 5);
    border_width(vbox, 5);
    gtk_container_add(GTK_CONTAINER(hbox), vbox);
    lbl = gtk_label_new(_("Mail setup"));
    gtk_notebook_append_page(GTK_NOTEBOOK(nb), hbox, lbl);  

    tbl = gtk_table_new(2, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_container_add(GTK_CONTAINER(vbox), tbl);

    lbl = gtk_label_new(_("SMTP server:"));
    gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
    gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, 0, 1, GTK_FILL, GTK_FILL, 0, 0);

    md.server_entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), md.server_entry, 1, 2, 0, 1);
    if (md.minfo->server != NULL) {
	gtk_entry_set_text(GTK_ENTRY(md.server_entry), md.minfo->server);
    }    

    lbl = gtk_label_new(_("port:"));
    gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
    gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, 1, 2, GTK_FILL, GTK_FILL, 0, 0);

    md.port_entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), md.port_entry, 1, 2, 1, 2);
    port_str = g_strdup_printf("%d", md.minfo->port);
    gtk_entry_set_text(GTK_ENTRY(md.port_entry), port_str);
    g_free(port_str);

    /* Create the "OK" button */
    md.ok = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(md.ok, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(md.dlg)->action_area), 
		       md.ok, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(md.ok), "clicked", 
		     G_CALLBACK(finalize_mail_settings), &md);
    gtk_widget_grab_default(md.ok);

    /* And a Cancel button */
    md.cancel = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(md.cancel, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(md.dlg)->action_area), 
		       md.cancel, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(md.cancel), "clicked", 
		     G_CALLBACK(finalize_mail_settings), &md);

    gtk_widget_set_size_request(md.dlg, 420, -1);
    gtk_widget_show_all(md.dlg);

    if (md.minfo->server == NULL) {
	gtk_notebook_set_current_page(GTK_NOTEBOOK(nb), 1);
    }

    gtk_window_set_modal(GTK_WINDOW(md.dlg), TRUE);
    gtk_main();

    return err;
}

#if GTK_MAJOR_VERSION < 2

static void mail_infobox (const char *msg, int err) 
{
    GtkWidget *w, *label, *button, *vbox, *hbox;

    w = gtk_window_new(GTK_WINDOW_DIALOG);

    gtk_container_border_width(GTK_CONTAINER(w), 5);
    gtk_window_position (GTK_WINDOW(w), GTK_WIN_POS_MOUSE);
    if (err) {
	gtk_window_set_title(GTK_WINDOW (w), _("gretl error"));
    } else {
	gtk_window_set_title(GTK_WINDOW (w), _("gretl info"));
    } 

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

#else /* GTK versions switch */

static void mail_infobox (const char *msg, int err)
{
    GtkWidget *dialog;

    dialog = gtk_message_dialog_new (NULL,
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     (err)? GTK_MESSAGE_ERROR : GTK_MESSAGE_INFO,
				     GTK_BUTTONS_CLOSE,
				     msg);
    gtk_dialog_run(GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
}

#endif

static int connect_to_server (char *hostname, unsigned short port) 
{
    gchar *msg;
    struct sockaddr_in soaddr;
    struct hostent *ip;
    int unit;

    ip = gethostbyname(hostname);
    if (ip == NULL) {
	msg = g_strdup_printf("Couldn't resolve name of server '%s': %s",
			      hostname, hstrerror(h_errno));
	mail_infobox(msg, 1);
	g_free(msg);
	return -1;
    }
    
    fprintf(stderr, "got server ip\n");

    unit = socket(PF_INET, SOCK_STREAM, 6); /* 6 = TCP */
    if (unit == -1) {
	mail_infobox("Couldn't open socket", 1);
	return -1;
    }

    fprintf(stderr, "got socket %d\n", unit);

    soaddr.sin_family = AF_INET;

    memcpy(&soaddr.sin_addr, &((struct in_addr*) ip->h_addr)->s_addr,
	   sizeof(struct in_addr));
    soaddr.sin_port = htons(port);

    if (connect(unit, (struct sockaddr*) &soaddr, sizeof soaddr) < 0) {
	msg = g_strdup_printf("Couldn't connect to %s", hostname);
	mail_infobox(msg, 1);
	g_free(msg);
	close(unit);
	return -1;
    } 

    fprintf(stderr, "connected socket OK\n");

    return unit;
}

static int
pop_before_smtp (const char *server, const char *sender)
{
    FILE *fp, *fq;
    char buf[1024];
    char pop_server[128];
    const char *p;
    int unit;

    p = strchr(server, '.');
    if (p == NULL) {
	return 1;
    }

    sprintf(pop_server, "pop%s", p);

    fprintf(stderr, "trying pop_before_smtp, with %s\n", pop_server);
    
    unit = connect_to_server(pop_server, 110);
    if (unit < 0) {
	return 1;
    } 

    fp = fdopen(unit, "r");
    fq = fdopen(unit, "w");

    if (fp == NULL || fq == NULL) {
	close(unit);
	return 1;
    }

    fgets(buf, sizeof buf, fp);
    fprintf(stderr, "buf: %s\n", buf);

    fputs("user cottrell\n", fq);
    fflush(fq);
    fgets(buf, sizeof buf, fp);
    fprintf(stderr, "buf: %s\n", buf);

    fputs("pass hozumi86\n", fq);   
    fflush(fq);
    fgets(buf, sizeof buf, fp);
    fprintf(stderr, "buf: %s\n", buf);
 
    fputs("quit\r\n", fq); 
    fflush(fq);
    fgets(buf, sizeof buf, fp);
    fprintf(stderr, "buf: %s\n", buf);
    
    fclose(fp);
    fclose(fq);
    
    close(unit);

    return 0;
}

static int
real_send_mail (FILE *infile, char *sender, char *recipient, 
		struct mail_info *minfo)
{
    char localhost[256] = "localhost";
    char buf[2048];
    gchar *errmsg;
    FILE *fp, *fq;
    int resp;
    int unit, err = 0;

    gethostname(localhost, sizeof localhost);

    unit = connect_to_server(minfo->server, minfo->port);
    if (unit < 0) {
	return 1;
    }

    fp = fdopen(unit, "r");
    fq = fdopen(unit, "w");

    if (fp == NULL || fq == NULL) {
	close(unit);
	return 1;
    }

    /* initial flush of stream */
    fgets(buf, sizeof buf, fp);

    fprintf(fq, "HELO %s\r\n", localhost);
    fflush(fq);

    fgets(buf, sizeof buf, fp);
    fprintf(stderr, "server response to HELO:\n%s\n", buf);

    fprintf(fq, "MAIL FROM: %s\r\n", sender);
    fflush(fq);

    fprintf(stderr, "sent MAIL FROM:\n", buf);

    fgets(buf, sizeof buf, fp);
    resp = atoi(buf);
    if (resp != 250) {
	chopstr(buf);
	errmsg = g_strdup_printf("Server response to MAIL FROM:\n%s", buf);
	mail_infobox(errmsg, 1);
	g_free(errmsg);
	err = 1;
	goto bailout;
    }

 try_again:

    fprintf(fq, "RCPT TO: %s\r\n", recipient);
    fflush(fq);

    fprintf(stderr, "sent RCPT TO:\n", buf);

    fgets(buf, sizeof buf, fp);
    resp = atoi(buf);
    if (resp == 553) {
	err = pop_before_smtp(minfo->server, sender);
	if (!err) {
	    goto try_again;
	}
    } else if (resp != 250 && resp != 251) {
	chopstr(buf);
	errmsg = g_strdup_printf("Server response to RCPT TO:\n%s", buf);
	mail_infobox(errmsg, 1);
	g_free(errmsg);
	err = 1;
	goto bailout;
    }

    fputs("DATA\r\n", fq);
    /* send composed message */
    while (fgets(buf, sizeof buf, infile)) {
	fputs(buf, fq);
    }
    fprintf(fq, "\n.\nQUIT\n");
    fflush(fq);

    fgets(buf, sizeof buf, fp);
    fprintf(stderr, "on QUIT, response: '%s'\n", buf);

 bailout:

    fclose(fp);
    fclose(fq);

    close(unit);

    return err;
}

static int pack_and_mail (const char *fname, struct msg_info *msg,
			  struct mail_info *minfo, char *tmpfname)
{
    const char *ctype;
    FILE *fp;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	perror(fname);
	err = 1;
    }

    if (is_data_file(fname)) {
	ctype = "application/x-gretldata";
    } else {
	ctype = "application/x-gretlscript";
    }

    if (!err) {
	err = encode(fp, fname, msg->note, msg->subj, msg->recip,
		     msg->sender, ctype, tmpfname);
    }

    if (!err) {
	fp = gretl_fopen(tmpfname, "r");
	if (fp == NULL) {
	    perror(tmpfname);
	    err = 1;
	} 
    }

    if (!err) {
	real_send_mail(fp, msg->sender, msg->recip, minfo);
	fclose(fp);
    }

    remove(tmpfname);

    return err;
}

int email_file (const char *fname, const char *userdir)
{
    struct mail_info minfo;
    struct msg_info msg;
    char temp[FILENAME_MAX];
    gchar *errmsg = NULL;
    int mval, err = 0;

    mail_info_init(&minfo);
    msg_init(&msg);

    sprintf(temp, "%smpack.XXXXXX", userdir);
    if (mktemp(temp) == NULL) {
	err = 1;
    }

    if (!err) {
	mval = mail_to_dialog(fname, &minfo, &msg);
	if (mval == MAIL_NO_RECIPIENT) {
	    errmsg = g_strdup("No address was given");
	} else if (mval == MAIL_NO_SERVER) {
	    errmsg = g_strdup("No SMTP was given");
	} else if (mval == MAIL_NO_SENDER) {
	    errmsg = g_strdup("No sender address was given");
	} else if (mval == MAIL_OK) {
	    err = pack_and_mail(fname, &msg, &minfo, temp);
	}
    }

    if (errmsg != NULL) {
	mail_infobox(errmsg, 1);
	g_free(errmsg);
	err = 1;
    }

    free_msg(&msg);
    free_mail_info(&minfo);

    return err;
}

