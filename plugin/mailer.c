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

/* Mailer plugin for gretl.  MIME packing is based on mpack, by John
   G. Myers.  Please see the files in ./mpack for the Carnegie Mellon
   copyright notice. */

#include "libgretl.h"
#include "version.h"

#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# include "gtk_compat.h"
#endif

#include "mpack/mpack.h"

typedef enum {
    MAIL_OK,
    MAIL_NO_RECIPIENT,
    MAIL_NO_SERVER,
    MAIL_NO_SENDER,
    MAIL_NO_PASS,
    MAIL_NO_MEM,
    MAIL_CANCEL
} MailError;

typedef enum {
    SMTP_OK,
    SMTP_NO_CONNECT,
    SMTP_NO_RELAY,
    SMTP_POP_FIRST,
    SMTP_BAD_SENDER,
    SMTP_BAD_ADDRESS,
    SMTP_OLD_SERVER,
    SMTP_ERR
} SMTPError;

typedef enum {
    SMTP_EHLO,
    SMTP_MAIL,
    SMTP_RCPT,
    SMTP_DATA,
    SMTP_DOT,
    SMTP_QUIT
} SMTPCode;

#define SBSIZE 4096

struct msg_info {
    char *recip;
    char *sender;
    char *subj;
    char *note;  
};  

struct mail_info {
    int err;
    char *sender;
    char *sig;
    int want_sig;
    char *server;
    unsigned short port;
    char *pop_server;
    char *pop_user;
    char *pop_pass;
    char *addrfile;
    GList *addrs;
};

struct mail_dialog {
    GtkWidget *dlg;
    GtkWidget *recip_combo;
    GtkWidget *reply_entry;
    GtkWidget *subj_entry;
    GtkWidget *note_entry;
    GtkWidget *server_entry;
    GtkWidget *port_entry;
    struct mail_info *minfo;
    struct msg_info *msg;
};

struct pop_dialog {
    GtkWidget *dlg;
    GtkWidget *server_entry;
    GtkWidget *user_entry;
    GtkWidget *pass_entry;
    struct mail_info *minfo;
};

static struct msg_info *mail_msg_new (void)
{
    struct msg_info *msg;

    msg = malloc(sizeof *msg);
    if (msg == NULL) {
	return NULL;
    }

    msg->recip = NULL;
    msg->sender = NULL;
    msg->subj = NULL;
    msg->note = NULL;

    return msg;
}

static void mail_msg_free (struct msg_info *msg)
{
    if (msg != NULL) {
	free(msg->recip);
	free(msg->sender);
	free(msg->subj);
	free(msg->note);
	free(msg);
    }
}

static struct mail_info *mail_info_new (void)
{
    struct mail_info *minfo;

    minfo = malloc(sizeof *minfo);
    if (minfo == NULL) {
	return NULL;
    }

    minfo->err = 0;
    minfo->sender = NULL;
    minfo->sig = NULL;
    minfo->want_sig = 1;
    minfo->server = NULL;
    minfo->port = 25;
    minfo->pop_server = NULL;
    minfo->pop_user = NULL;
    minfo->pop_pass = NULL;
    minfo->addrfile = NULL;
    minfo->addrs = NULL;

    return minfo;
}

static void mail_info_free (struct mail_info *minfo)
{
    GList *tmp;

    if (minfo == NULL) {
	return;
    }

    free(minfo->sender);
    free(minfo->sig);
    free(minfo->server);
    free(minfo->pop_server);
    free(minfo->pop_user);
    free(minfo->pop_pass);
    free(minfo->addrfile);

    tmp = minfo->addrs;
    while (tmp != NULL) {
	g_free(tmp->data);
	tmp = g_list_next(tmp);
    }

    free(minfo);
}

static void save_email_info (struct mail_info *minfo)
{
    FILE *fp;

    fp = gretl_fopen(minfo->addrfile, "w");

    if (fp != NULL) {
	GList *list = minfo->addrs;
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
	if (minfo->pop_server != NULL && *minfo->pop_server != '\0') {
	    fprintf(fp, "POP server: %s\n", minfo->pop_server);
	}
	if (minfo->pop_user != NULL && *minfo->pop_user != '\0') {
	    fprintf(fp, "POP user: %s\n", minfo->pop_user);
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

static void cancel_mail_settings (GtkWidget *w, struct mail_dialog *md)
{
    md->minfo->err = MAIL_CANCEL;
    gtk_widget_destroy(md->dlg);
}

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 6)

static gchar *my_combo_box_get_active_text (GtkComboBox *box)
{
    GtkWidget *w = gtk_bin_get_child(GTK_BIN(box));
    gchar *ret = NULL;

    if (GTK_IS_ENTRY(w)) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));

	if (s != NULL) {
	    ret = g_strdup(s);
	}
    } 

    return ret;
}

#endif

static void finalize_mail_settings (GtkWidget *w, struct mail_dialog *md)
{
    GList *list = NULL;
    struct mail_info *minfo = md->minfo;
    struct msg_info *msg = md->msg;
    gchar *recip;
    const gchar *txt;
    int err = MAIL_OK;
    int save = 0;

    list = minfo->addrs;

    /* recipient */
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 6)
    recip = my_combo_box_get_active_text(GTK_COMBO_BOX(md->recip_combo));
#elif GTK_MAJOR_VERSION < 3
    recip = gtk_combo_box_get_active_text(GTK_COMBO_BOX(md->recip_combo));
#else
    recip = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(md->recip_combo));    
#endif

    if (recip != NULL && *recip != '\0') {
	int i = 0;

	msg->recip = g_strdup(recip);
	fprintf(stderr, "targ = '%s'\n", msg->recip);
	save = 1;
	while (list) {
	    if (!strcmp(recip, (char *) list->data)) {
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
	    minfo->addrs = g_list_prepend(minfo->addrs, g_strdup(recip));
	} 
    } else {
	err = MAIL_NO_RECIPIENT;
    }

    g_free(recip);

    if (!err) {
	/* reply-to address */
	txt = gtk_entry_get_text(GTK_ENTRY(md->reply_entry));
	if (txt != NULL && *txt != '\0') {
	    msg->sender = g_strdup(txt);
	    if (minfo->sender == NULL) {
		minfo->sender = g_strdup(txt);
		save = 1;
	    } else if (strcmp(txt, minfo->sender)) {
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

    md->minfo->err = err;

    if (save) {
	save_email_info(minfo);
    }

    gtk_widget_destroy(md->dlg);
}

static void cancel_pop_settings (GtkWidget *w, struct pop_dialog *pd)
{
    pd->minfo->err = MAIL_CANCEL;
    gtk_widget_destroy(pd->dlg);
}

static void finalize_pop_settings (GtkWidget *w, struct pop_dialog *pd)
{
    struct mail_info *minfo = pd->minfo;
    const gchar *txt;
    int err = MAIL_OK;

    if (!err) {
	/* server */
	txt = gtk_entry_get_text(GTK_ENTRY(pd->server_entry));
	if (txt != NULL && *txt != '\0') {
	    minfo->pop_server = g_strdup(txt);
	    fprintf(stderr, "POP server = '%s'\n", minfo->pop_server);
	} else {
	    err = MAIL_NO_SERVER;
	}
    }

    if (!err) {
	/* username */
	txt = gtk_entry_get_text(GTK_ENTRY(pd->user_entry));
	if (txt != NULL && *txt != '\0') {
	    minfo->pop_user = g_strdup(txt);
	    fprintf(stderr, "username = '%s'\n", minfo->pop_user);
	} else {
	    err = MAIL_NO_SENDER;
	}
    } 

    if (!err) {
	/* password */
	txt = gtk_entry_get_text(GTK_ENTRY(pd->pass_entry));
	if (txt != NULL && *txt != '\0') {
	    minfo->pop_pass = g_strdup(txt);
	    fprintf(stderr, "got %d character password\n", (int) strlen(txt));
	} else {
	    err = MAIL_NO_PASS;
	}
    }

    if (!err) {
	save_email_info(minfo);
    }

    pd->minfo->err = err;

    gtk_widget_destroy(pd->dlg);
}


static void border_width (GtkWidget *w, int b)
{
    gtk_container_set_border_width(GTK_CONTAINER(w), b); 
}

static void set_dialog_border_widths (GtkWidget *dlg)
{
    GtkWidget *box;
    int w1 = 10, w2 = 5;

    box = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    border_width(box, w1);
    gtk_box_set_spacing(GTK_BOX(box), w2);
    box = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    border_width(box, w2);
}

static void get_email_info (struct mail_info *minfo)
{
    GList *addrs = NULL;
    FILE *fp;

    minfo->addrfile = g_strdup_printf("%sgretl.addresses", gretl_dotdir());

    fp = gretl_fopen(minfo->addrfile, "r");
    if (fp != NULL) {
	char line[128];

	while (fgets(line, sizeof line, fp)) {
	    if (string_is_blank(line)) {
		continue;
	    }
	    gretl_strstrip(line);
	    if (!strncmp(line, "Reply-To:", 9)) {
		minfo->sender = g_strdup(line + 10);
	    } else if (!strncmp(line, "SMTP server:", 12)) {
		minfo->server = g_strdup(line + 13);
	    } else if (!strncmp(line, "SMTP port:", 10)) {
		minfo->port = atoi(line + 11);
	    } else if (!strncmp(line, "POP server:", 11)) {
		minfo->pop_server = g_strdup(line + 12);
	    } else if (!strncmp(line, "POP user:", 9)) {
		minfo->pop_user = g_strdup(line + 10);
	    } else {
		addrs = g_list_append(addrs, g_strdup(line));
	    }
	}

	fclose(fp);
    } 

    minfo->addrs = addrs;
}

static gboolean 
mail_dialog_delete (GtkWidget *w, GdkEvent *e, struct mail_info *minfo)
{
    minfo->err = MAIL_CANCEL;
    return FALSE;
}

static int is_data_file (const char *fname)
{
    int ret = 1;

    if (fname != NULL && strlen(fname) > 4) {
	const char *ext = strrchr(fname, '.');

	if (ext != NULL) {
	    if (!strcmp(ext, ".inp")) {
		ret = 0;
	    } else if (!strcmp(ext, ".gfn")) {
		ret = 0;
	    } else if (!strcmp(ext, ".gretl")) {
		ret = 0;
	    }
	}	    
    }

    return ret;
}

static void sig_callback (GtkWidget *w, struct mail_info *minfo)
{
    minfo->want_sig = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
}

static void mail_dialog_quit (GtkWidget *w, gpointer p)
{
    gtk_main_quit();
}

static void set_combo_strings_from_list (GtkComboBox *box, GList *list)
{
    GList *mylist = list;

    while (mylist != NULL) {
#if GTK_MAJOR_VERSION >= 3
	gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(box), 
				       mylist->data);
#else
	gtk_combo_box_append_text(box, mylist->data);
#endif
	mylist = mylist->next;
    }
}

static void 
mail_to_dialog (const char *fname, struct mail_info *minfo, 
		struct msg_info *msg)
{
    const gchar *lbls[] = {
	N_("To:"),
	N_("Reply-To:"),
	N_("Subject:"),
	N_("Note:")
    };
    GtkWidget *tbl, *lbl, *vbox;
    GtkWidget *nb, *hbox, *button;
    gchar *port_str;
    const char *short_fname, *p;
    struct mail_dialog md;
    int i, datafile, nrows;

    md.dlg = gtk_dialog_new();
    md.minfo = minfo;
    md.msg = msg;

    get_email_info(minfo);
    minfo->sig = get_signature();
    minfo->want_sig = minfo->sig != NULL;

    g_signal_connect(G_OBJECT(md.dlg), "delete_event", 
		     G_CALLBACK(mail_dialog_delete), minfo);

    g_signal_connect(G_OBJECT(md.dlg), "destroy", 
		     G_CALLBACK(mail_dialog_quit), NULL);

    gtk_window_set_title(GTK_WINDOW(md.dlg), _("gretl: send mail"));
    set_dialog_border_widths(md.dlg);
#if GTK_MAJOR_VERSION < 3    
    gtk_dialog_set_has_separator(GTK_DIALOG(md.dlg), FALSE);
#endif    
    gtk_window_set_position(GTK_WINDOW(md.dlg), GTK_WIN_POS_MOUSE);

    nb = gtk_notebook_new();
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(md.dlg));
    gtk_container_add(GTK_CONTAINER(vbox), nb);
    hbox = gtk_hbox_new(FALSE, 5);
    border_width(hbox, 5);
    vbox = gtk_vbox_new(FALSE, 5);
    border_width(vbox, 5);
    gtk_container_add(GTK_CONTAINER(hbox), vbox);
    lbl = gtk_label_new(_("Message"));
    gtk_notebook_append_page(GTK_NOTEBOOK(nb), hbox, lbl);    

    nrows = (minfo->sig == NULL)? 4 : 5;

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
#if GTK_MAJOR_VERSION >= 3
	    w = gtk_combo_box_text_new_with_entry();
#else
	    w = gtk_combo_box_entry_new_text();
#endif
	    if (minfo->addrs != NULL) {
		set_combo_strings_from_list(GTK_COMBO_BOX(w), minfo->addrs);
	    } 
	} else {
	    w = gtk_entry_new();
	}

	if (i == 1) {
	    if (minfo->sender != NULL) {
		gtk_entry_set_text(GTK_ENTRY(w), minfo->sender);
	    }
	} else if (i == 2) {
	    gtk_entry_set_text(GTK_ENTRY(w), (datafile)? "dataset" : "script");
	} else if (i == 3) {
	    gchar *note;

	    if (datafile) {
		note = g_strdup_printf(_("Please find the gretl data file %s attached."),
				       short_fname);
	    } else {
		note = g_strdup_printf(_("Please find the gretl script %s attached."),
				       short_fname);
	    }		
	    gtk_entry_set_text(GTK_ENTRY(w), note);
	    g_free(note);
	}

	if (i == 0) {
	    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(w));

	    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
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

    if (minfo->sig != NULL) {
	GtkWidget *w;

	w = gtk_check_button_new_with_label(_("Append signature"));
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

    tbl = gtk_table_new(2, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_container_add(GTK_CONTAINER(vbox), tbl);

    lbl = gtk_label_new(_("SMTP server:"));
    gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
    gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, 0, 1, GTK_FILL, GTK_FILL, 0, 0);

    md.server_entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), md.server_entry, 1, 3, 0, 1);
    if (minfo->server != NULL) {
	gtk_entry_set_text(GTK_ENTRY(md.server_entry), minfo->server);
    }    

    lbl = gtk_label_new(_("port:"));
    gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
    gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, 1, 2, GTK_FILL, GTK_FILL, 0, 0);

    md.port_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(md.port_entry), 5);
    gtk_entry_set_width_chars(GTK_ENTRY(md.port_entry), 8);

    gtk_table_attach_defaults(GTK_TABLE(tbl), md.port_entry, 1, 2, 1, 2);
    port_str = g_strdup_printf("%d", minfo->port);
    gtk_entry_set_text(GTK_ENTRY(md.port_entry), port_str);
    g_free(port_str);
    lbl = gtk_label_new("                     ");
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 2, 3, 1, 2);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(md.dlg));

    /* Cancel button */
    button = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(cancel_mail_settings), &md);

    /* Create the "OK" button */
    button = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(finalize_mail_settings), &md);
    gtk_widget_grab_default(button);

    gtk_widget_set_size_request(md.dlg, 420, -1);
    gtk_widget_show_all(md.dlg);

    if (minfo->server == NULL) {
	gtk_notebook_set_current_page(GTK_NOTEBOOK(nb), 1);
    }

    gtk_window_set_modal(GTK_WINDOW(md.dlg), TRUE);
    gtk_main();
}

static void pop_info_dialog (struct mail_info *minfo)
{
    const gchar *lbls[] = {
	N_("POP server:"),
	N_("Username:"),
	N_("Password:")
    };
    GtkWidget *tbl, *lbl, *button;
    GtkWidget *hbox, *vbox;
    struct pop_dialog pd;
    int i;

    pd.dlg = gtk_dialog_new();
    pd.minfo = minfo;

    g_signal_connect(G_OBJECT(pd.dlg), "delete_event", 
		     G_CALLBACK(mail_dialog_delete), minfo);

    g_signal_connect(G_OBJECT(pd.dlg), "destroy", 
		     G_CALLBACK(mail_dialog_quit), NULL);

    gtk_window_set_title(GTK_WINDOW(pd.dlg), _("gretl: POP info"));
    set_dialog_border_widths(pd.dlg);
    gtk_window_set_position(GTK_WINDOW(pd.dlg), GTK_WIN_POS_MOUSE);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(pd.dlg));

    tbl = gtk_table_new(3, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_container_add(GTK_CONTAINER(vbox), tbl);

    for (i=0; i<3; i++) {
	GtkWidget *w;

	lbl = gtk_label_new(_(lbls[i]));
	gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
	gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, i, i+1, GTK_FILL, GTK_FILL, 0, 0);

	w = gtk_entry_new();

	if (i == 0) {
	    if (minfo->pop_server != NULL) {
		gtk_entry_set_text(GTK_ENTRY(w), minfo->pop_server);
	    }
	} else if (i == 1) {
	    if (minfo->pop_user != NULL) {
		gtk_entry_set_text(GTK_ENTRY(w), minfo->pop_user);
	    }
	} else if (i == 2) {
	    if (minfo->pop_pass != NULL) {
		gtk_entry_set_text(GTK_ENTRY(w), minfo->pop_pass);
	    }
	    gtk_entry_set_visibility(GTK_ENTRY(w), FALSE);
	}

	gtk_entry_set_activates_default(GTK_ENTRY(w), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), w, 1, 2, i, i+1);

	if (i == 0) {
	    pd.server_entry = w;
	} else if (i == 1) {
	    pd.user_entry = w;
	} else if (i == 2) {
	    pd.pass_entry = w;
	} 
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(pd.dlg));

    /* Cancel button */
    button = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(cancel_pop_settings), &pd);

    /* "OK" button */
    button = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(finalize_pop_settings), &pd);
    gtk_widget_grab_default(button);

    gtk_widget_set_size_request(pd.dlg, 360, -1);
    gtk_widget_show_all(pd.dlg);

    gtk_window_set_modal(GTK_WINDOW(pd.dlg), TRUE);
    gtk_main();
}

static void mail_infobox (const char *msg, int err)
{
    GtkWidget *dialog;

    dialog = gtk_message_dialog_new (NULL,
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     (err)? GTK_MESSAGE_ERROR : GTK_MESSAGE_INFO,
				     GTK_BUTTONS_CLOSE,
				     "%s",
				     msg);
    gtk_dialog_run(GTK_DIALOG (dialog));
    gtk_widget_destroy(dialog);
}

#define VERBOSE 1

static int get_server_response (int fd, char *buf)
{
    int ret;

    memset(buf, 0, SBSIZE);
#if VERBOSE
    fputs("doing read() on socket...\n", stderr);
#endif
    ret = read(fd, buf, SBSIZE - 1);
#if VERBOSE
    fprintf(stderr, "response:\n%s\n", buf);
#endif
    return ret;
}

static int send_to_server (FILE *fp, const char *template, ...)
{
    va_list args;
    int plen = 0;

#if VERBOSE
    char word[32] = {0};

    sscanf(template, "%31s", word);
    fprintf(stderr, "sending %s...\n", word);
#endif

    va_start(args, template);
    plen = vfprintf(fp, template, args);
    va_end(args);

    fflush(fp);

    return plen;
}

#ifndef HAVE_IN_ADDR
struct in_addr {
    unsigned long s_addr;
}; 
#endif

#ifndef HAVE_SOCKADDR_IN
struct sockaddr_in {
    short int          sin_family;
    unsigned short int sin_port;
    struct in_addr     sin_addr;
    unsigned char      sin_zero[8];
};
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

    soaddr.sin_family = AF_INET;
    memcpy(&soaddr.sin_addr, &((struct in_addr *) ip->h_addr)->s_addr,
	   sizeof(struct in_addr));
    soaddr.sin_port = htons(port);
    memset(&soaddr.sin_zero, '\0', 8);

    if (connect(unit, (struct sockaddr *) &soaddr, sizeof soaddr) < 0) {
	msg = g_strdup_printf("Couldn't connect to %s", hostname);
	mail_infobox(msg, 1);
	g_free(msg);
	close(unit);
	return -1;
    } 

    return unit;
}

static int set_pop_defaults (struct mail_info *minfo)
{
    char *p;

    if (minfo->server == NULL || minfo->sender == NULL) {
	/* these must be defined at this point */
	return 1;
    }

    if (minfo->pop_server == NULL) {
	p = strchr(minfo->server, '.');
	if (p != NULL) {
	    minfo->pop_server = g_strdup_printf("pop%s", p);
	}
    }

    if (minfo->pop_user == NULL) {
	p = strchr(minfo->sender, '@');
	if (p != NULL) {
	    minfo->pop_user = g_strdup(minfo->sender);
	    p = strchr(minfo->pop_user, '@');
	    *p = '\0';
	}
    }

    return 0;
}

static int get_POP_error (char *buf)
{
    int err = 0;

    if (*buf == '-') {
	gchar *errmsg;

	gretl_strstrip(buf);
	errmsg = g_strdup_printf("POP server said:\n%s", buf);
	mail_infobox(errmsg, 1);
	g_free(errmsg);
	err = 1;
    }

    return err;
}

static int pop_login (struct mail_info *minfo)
{
    FILE *fp;
    char buf[SBSIZE];
    int unit, err;

    set_pop_defaults(minfo);
    pop_info_dialog(minfo);

    if (minfo->err) {
	return 1;
    }

    fprintf(stderr, "trying POP before SMTP, with %s\n", minfo->pop_server);
    
    unit = connect_to_server(minfo->pop_server, 110);
    if (unit < 0) {
	return 1;
    } 

    fp = fdopen(unit, "w");
    if (fp == NULL) {
	close(unit);
	return 1;
    }

    get_server_response(unit, buf);

    send_to_server(fp, "USER %s\n", minfo->pop_user);
    get_server_response(unit, buf);
    err = get_POP_error(buf);

    if (!err) {
	send_to_server(fp, "PASS %s\n", minfo->pop_pass);   
	get_server_response(unit, buf);
	err = get_POP_error(buf);
    }
 
    send_to_server(fp, "QUIT\r\n"); 
    get_server_response(unit, buf);
    
    fclose(fp);
    close(unit);

    return err;
}

static int get_SMTP_error (char *buf, SMTPCode code)
{
    gchar *errmsg = NULL;
    int resp = atoi(buf);
    int err = SMTP_OK;

    if (code == SMTP_EHLO) {
	if (resp == 500) {
	    err = SMTP_OLD_SERVER;
	} else if (resp != 250) {
	    gretl_strstrip(buf);
	    errmsg = g_strdup_printf("Server response to . :\n%s", buf);
	    err = SMTP_ERR;
	}
    } else if (code == SMTP_MAIL || code == SMTP_RCPT) {
	if (resp == 553 && strstr(buf, "must check")) {
	    err = SMTP_POP_FIRST;
	} else if (resp != 250) {
	    gretl_strstrip(buf);
	    errmsg = g_strdup_printf("Server response to RCPT:\n%s", buf);
	    err = SMTP_ERR;
	}
    } else if (code == SMTP_DATA) {
	if (resp != 354) {
	    gretl_strstrip(buf);
	    errmsg = g_strdup_printf("Server response to RCPT:\n%s", buf);
	    err = SMTP_ERR;
	}
    } else if (code == SMTP_DOT) {
	if (resp != 250) {
	    gretl_strstrip(buf);
	    errmsg = g_strdup_printf("Server response to . :\n%s", buf);
	    err = SMTP_ERR;
	}
    } 

    if (errmsg != NULL) {
	mail_infobox(errmsg, 1);
	g_free(errmsg);
    }

    return err;
}

static int
smtp_send_mail (FILE *infile, char *sender, char *recipient, 
		struct mail_info *minfo)
{
    char localhost[256] = "localhost";
    char buf[SBSIZE];
    FILE *fp;
    int unit, err = SMTP_OK;

    gethostname(localhost, sizeof localhost);
    fprintf(stderr, "localhost = '%s'\n", localhost);

    unit = connect_to_server(minfo->server, minfo->port);
    if (unit < 0) {
	return SMTP_NO_CONNECT;
    }

    fprintf(stderr, "opened SMTP socket, unit = %d\n", unit);

    fp = fdopen(unit, "w");
    if (fp == NULL) {
	close(unit);
	return SMTP_ERR;
    }

    get_server_response(unit, buf);

    send_to_server(fp, "EHLO %s\r\n", localhost);
    get_server_response(unit, buf);
    err = get_SMTP_error(buf, SMTP_EHLO);
    if (err == SMTP_OLD_SERVER) {
	send_to_server(fp, "HELO %s\r\n", localhost);
	get_server_response(unit, buf);
	err = get_SMTP_error(buf, SMTP_EHLO);
    }
    if (err) goto bailout;

    send_to_server(fp, "MAIL FROM:<%s>\r\n", sender);
    get_server_response(unit, buf);
    err = get_SMTP_error(buf, SMTP_MAIL);
    if (err) goto bailout;

    send_to_server(fp, "RCPT TO:<%s>\r\n", recipient);
    get_server_response(unit, buf);
    err = get_SMTP_error(buf, SMTP_RCPT);
    if (err) goto bailout;

    send_to_server(fp, "DATA\r\n");
    get_server_response(unit, buf);
    err = get_SMTP_error(buf, SMTP_DATA);
    if (err) goto bailout;

#if VERBOSE
    fputs("sending actual message...\n", stderr);
#endif    
    while (fgets(buf, sizeof buf - 1, infile)) {
	int n = strlen(buf);

	/* rfc2821: ensure CRLF termination */
	if (buf[n-1] == '\n' && buf[n-2] != '\r') {
	    buf[n-1] = '\r';
	    buf[n] = '\n';
	    buf[n+1] = '\0';
	}
	fputs(buf, fp);
    }
    fputs("\r\n.\r\n", fp);
    fflush(fp);

    get_server_response(unit, buf);
    err = get_SMTP_error(buf, SMTP_DOT);

 bailout:

    send_to_server(fp, "QUIT\r\n");
    get_server_response(unit, buf);

    fclose(fp);
    close(unit);

    return err;
}

static int pack_and_mail (const char *fname, 
			  struct msg_info *msg,
			  struct mail_info *minfo, 
			  const char *userdir)
{
    char tmpfname[FILENAME_MAX];
    FILE *fpin, *fpout;
    int err = 0;

    fpin = gretl_fopen(fname, "rb");
    if (fpin == NULL) {
	perror(fname);
	err = 1;
    }

    sprintf(tmpfname, "%smpack.XXXXXX", userdir);
    fpout = gretl_mktemp(tmpfname, "wb");
    if (fpout == NULL) {
	err = 1;
    }

    if (!err) {
	const char *ctype = is_data_file(fname) ? 
	    "application/x-gretldata" : "application/x-gretlscript";

	err = mpack_encode(fpin, fname, msg->note, msg->subj, msg->recip,
			   msg->sender, ctype, fpout);
    }

    if (fpin != NULL) {
	fclose(fpin);
    }

    if (fpout != NULL) {
	fclose(fpout);
    }

    if (!err) {
	fpin = gretl_fopen(tmpfname, "rb");
	if (fpin == NULL) {
	    perror(tmpfname);
	    err = 1;
	} else {
	    err = smtp_send_mail(fpin, msg->sender, msg->recip, minfo);
	    if (err == SMTP_POP_FIRST) {
		err = pop_login(minfo);
		if (!err) {
		    err = smtp_send_mail(fpin, msg->sender, msg->recip, minfo);
		}
	    }
	    fclose(fpin);
	}
    }

    gretl_remove(tmpfname);

    return err;
}

int email_file (const char *fname, const char *userdir)
{
    struct mail_info *minfo = NULL;
    struct msg_info *msg = NULL;
    gchar *errmsg = NULL;
    int err = 0;

    minfo = mail_info_new();
    msg = mail_msg_new();

    if (minfo == NULL || msg == NULL) {
	mail_msg_free(msg);
	mail_info_free(minfo);
	return E_ALLOC;
    }

    mail_to_dialog(fname, minfo, msg);

    if (minfo->err == MAIL_NO_RECIPIENT) {
	errmsg = g_strdup("No address was given");
    } else if (minfo->err == MAIL_NO_SERVER) {
	errmsg = g_strdup("No SMTP was given");
    } else if (minfo->err == MAIL_NO_SENDER) {
	errmsg = g_strdup("No sender address was given");
    } else if (minfo->err == MAIL_NO_MEM) {
	errmsg = g_strdup("Out of memory");
    } else if (minfo->err == MAIL_OK) {
	err = pack_and_mail(fname, msg, minfo, userdir);
    }

    if (errmsg != NULL) {
	mail_infobox(errmsg, 1);
	g_free(errmsg);
	err = 1;
    }

    mail_msg_free(msg);
    mail_info_free(minfo);

    return err;
}

