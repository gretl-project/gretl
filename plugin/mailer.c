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
#include "gretl_www.h"

#include "../pixmaps/eye.xpm"

#include <gtk/gtk.h>
#include <stdlib.h>

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# include "gtk_compat.h"
#endif

#include "mpack/mpack.h"

struct msg_info {
    gchar *recip;
    gchar *sender;
    gchar *subj;
    gchar *note;
    int canceled;
    int changed;
};

struct mail_info {
    gchar *sender;
    gchar *sig;
    int want_sig;
    gchar *server;
    gchar *mail_user;
    gchar *mail_pass;
    int pass_needed;
    gchar *addrfile;
    GList *addrs;
    int store;
};

struct mail_dialog {
    GtkWidget *dlg;
    GtkWidget *nb;
    GtkWidget *recip_combo;
    GtkWidget *sender_entry;
    GtkWidget *subj_entry;
    GtkWidget *note_entry;
    GtkWidget *server_entry;
    GtkWidget *user_entry;
    GtkWidget *pass_check;
    struct mail_info *minfo;
    struct msg_info *msg;
};

struct passwd_dialog {
    GtkWidget *dlg;
    GtkWidget *entry;
    struct mail_info *minfo;
    struct msg_info *msg;
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
    msg->canceled = 0;
    msg->changed = 0;

    return msg;
}

static void mail_msg_free (struct msg_info *msg)
{
    if (msg != NULL) {
	g_free(msg->recip);
	g_free(msg->sender);
	g_free(msg->subj);
	g_free(msg->note);
	g_free(msg);
    }
}

static struct mail_info *mail_info_new (void)
{
    struct mail_info *minfo;

    minfo = malloc(sizeof *minfo);
    if (minfo == NULL) {
	return NULL;
    }

    minfo->sender = NULL;
    minfo->sig = NULL;
    minfo->want_sig = 1;
    minfo->server = NULL;
    minfo->mail_user = NULL;
    minfo->mail_pass = NULL;
    minfo->pass_needed = 1;
    minfo->addrfile = NULL;
    minfo->addrs = NULL;
    minfo->store = 0;

    return minfo;
}

static void mail_info_free (struct mail_info *minfo)
{
    GList *tmp;

    if (minfo == NULL) {
	return;
    }

    g_free(minfo->sender);
    g_free(minfo->sig);
    g_free(minfo->server);
    g_free(minfo->mail_user);
    g_free(minfo->mail_pass);
    g_free(minfo->addrfile);

    tmp = minfo->addrs;
    while (tmp != NULL) {
	g_free(tmp->data);
	tmp = g_list_next(tmp);
    }

    free(minfo);
}

static void dump_content (gchar *s, FILE *fp)
{
    guchar u[2];
    gint n = strlen(s) + 1;
    int i, k = 13;

    fwrite(&n, sizeof n, 1, fp);
    for (i=0; i<n; i++) {
	u[0] = s[i] / k;
	u[1] = s[i] % k;
	fwrite(u, 1, 2, fp);
    }
}

static gchar *grab_content (FILE *fp)
{
    gchar *content = NULL;
    gint n, got;

    got = fread(&n, sizeof n, 1, fp);

    if (got == 1) {
	guchar u[2];
	int i, k = 13;

	content = g_malloc0(n);
	for (i=0; i<n; i++) {
	    got = fread(u, 1, 2, fp);
	    if (got == 2) {
		content[i] = k * u[0] + u[1];
	    }
	}
    }

    return content;
}

static void write_mail_info (struct mail_info *minfo)
{
    FILE *fp;

    fp = gretl_fopen(minfo->addrfile, "wb");

    if (fp != NULL) {
	GString *str = g_string_new(NULL);
	GList *list = minfo->addrs;
	gchar *content = NULL;
	int i, maxaddrs = 10;

	if (minfo->sender != NULL && *minfo->sender != '\0') {
	    g_string_append_printf(str, "Reply-To: %s\n", minfo->sender);
	}
	if (minfo->server != NULL && *minfo->server != '\0') {
	    g_string_append_printf(str, "SMTP server: %s\n", minfo->server);
	}
	if (minfo->mail_user != NULL && *minfo->mail_user != '\0') {
	    g_string_append_printf(str, "mail user: %s\n", minfo->mail_user);
	}
	if (minfo->store && minfo->mail_pass != NULL && *minfo->mail_pass != '\0') {
	    g_string_append_printf(str, "mail pass: %s\n", minfo->mail_pass);
	}
	g_string_append_printf(str, "password needed: %d\n", minfo->pass_needed);
	for (i=0; i<maxaddrs && list != NULL; i++) {
	    g_string_append_printf(str, "%s\n", (char *) list->data);
	    list = list->next;
	}

	content = g_string_free(str, FALSE);
	dump_content(content, fp);
	g_free(content);
	fclose(fp);
    }
}

static void read_mail_info (struct mail_info *minfo)
{
    gchar *content = NULL;
    GList *addrs = NULL;
    FILE *fp;

    minfo->addrfile = g_strdup_printf("%smail.dat", gretl_dotdir());
    fp = gretl_fopen(minfo->addrfile, "rb");

    if (fp != NULL) {
	content = grab_content(fp);
	fclose(fp);
    }

    if (content != NULL) {
	char line[128];

	bufgets_init(content);
	while (bufgets(line, sizeof line, content)) {
	    if (string_is_blank(line)) {
		continue;
	    }
	    gretl_strstrip(line);
	    if (!strncmp(line, "Reply-To:", 9)) {
		minfo->sender = g_strdup(line + 10);
	    } else if (!strncmp(line, "SMTP server:", 12)) {
		minfo->server = g_strdup(line + 13);
	    } else if (!strncmp(line, "mail user:", 10)) {
		minfo->mail_user = g_strdup(line + 11);
	    } else if (!strncmp(line, "mail pass:", 10)) {
		minfo->mail_pass = g_strdup(line + 11);
	    } else if (!strncmp(line, "password needed:", 16)) {
		minfo->pass_needed = atoi(line + 17);
	    } else {
		addrs = g_list_append(addrs, g_strdup(line));
	    }
	}
	bufgets_finalize(content);
	g_free(content);
    }

    /* default to google's SMTP server */
    if (minfo->server == NULL) {
	minfo->server = g_strdup("smtps://smtp.gmail.com:465");
    } else if (strncmp(minfo->server, "smtp://", 7) &&
	       strncmp(minfo->server, "smtps://", 8)) {
	g_free(minfo->server);
	minfo->server = g_strdup("smtps://smtp.gmail.com:465");
    }

    if (minfo->sender == NULL && minfo->mail_user != NULL &&
	strchr(minfo->mail_user, '@') != NULL) {
	minfo->sender = g_strdup(minfo->mail_user);
    }

    if (minfo->mail_pass != NULL) {
	minfo->store = 1;
    }

    minfo->addrs = addrs;
}

static char *get_signature (void)
{
    char *home = getenv("HOME");
    gchar *sig = NULL;

    if (home != NULL) {
	gchar *sigfile = g_strdup_printf("%s/.signature", home);
	FILE *fp;

	fp = gretl_fopen(sigfile, "r");
	if (fp != NULL) {
	    GString *tmp = g_string_new(NULL);
	    char line[128];

	    while (fgets(line, sizeof line, fp)) {
		tmp = g_string_append(tmp, line);
	    }
	    fclose(fp);
	    sig = g_string_free(tmp, 0);
	}
	g_free(sigfile);
    }

    return sig;
}

static int finalize_mail_settings (struct mail_dialog *md)
{
    struct mail_info *minfo = md->minfo;
    struct msg_info *msg = md->msg;
    GList *list = minfo->addrs;
    GtkWidget *targ = NULL;
    gchar *recip = NULL;
    const gchar *sender = NULL;
    const gchar *server = NULL;
    const gchar *user = NULL;
    const gchar *txt;
    int i, save = 0;

#if GTK_MAJOR_VERSION < 3
    recip = gtk_combo_box_get_active_text(GTK_COMBO_BOX(md->recip_combo));
#else
    recip = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(md->recip_combo));
#endif

    /* check for missing info */
    if (recip == NULL || *recip == '\0' || !strchr(recip, '@')) {
	targ = md->recip_combo;
    }
    if (targ == NULL) {
	sender = gtk_entry_get_text(GTK_ENTRY(md->sender_entry));
	if (sender == NULL || *sender == '\0' || !strchr(sender, '@')) {
	    targ = md->sender_entry;
	}
    }
    if (targ == NULL) {
	server = gtk_entry_get_text(GTK_ENTRY(md->server_entry));
	if (server == NULL || *server == '\0') {
	    targ = md->server_entry;
	}
    }
    if (targ == NULL) {
	user = gtk_entry_get_text(GTK_ENTRY(md->user_entry));
	if (user == NULL || *user == '\0') {
	    targ = md->user_entry;
	}
    }
    if (targ != NULL) {
	/* focus the widget holding the missing field */
	int page = gtk_notebook_page_num(GTK_NOTEBOOK(md->nb), targ);

	gtk_notebook_set_current_page(GTK_NOTEBOOK(md->nb), page);
	gtk_widget_grab_focus(targ);
	g_free(recip);
	return 1;
    }

    /* handle recipient info */
    msg->recip = recip;
    fprintf(stderr, "targ = '%s'\n", msg->recip);
    save = 1;
    i = 0;
    while (list != NULL) {
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

    /* sender's address */
    msg->sender = g_strdup(sender);
    if (minfo->sender == NULL || strcmp(sender, minfo->sender)) {
	g_free(minfo->sender);
	minfo->sender = g_strdup(sender);
	save = 1;
    }
    fprintf(stderr, "sender = '%s'\n", msg->sender);

    /* SMTP server */
    if (minfo->server == NULL || strcmp(server, minfo->server)) {
	g_free(minfo->server);
	minfo->server = g_strdup(server);
	save = 1;
    }
    fprintf(stderr, "server = '%s'\n", minfo->server);

    /* email userID */
    if (minfo->mail_user == NULL || strcmp(user, minfo->mail_user)) {
	g_free(minfo->mail_user);
	minfo->mail_user = g_strdup(user);
	save = 1;
    }
    fprintf(stderr, "username = '%s'\n", minfo->mail_user);

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
	    g_free(minfo->sig);
	    minfo->sig = NULL;
	}
	if (minfo->sig != NULL) {
	    msg->note = g_strdup_printf("%s\n--\n%s\n", txt, minfo->sig);
	} else {
	    msg->note = g_strdup_printf("%s\n", txt);
	}
    }

    if (save) {
	write_mail_info(minfo);
    }

    return 0;
}

/* callback for "response" signal from @dlg */

static void mail_dialog_callback (GtkDialog *dlg, gint id, int *ret)
{
    int err = 0;

    if (id == GTK_RESPONSE_REJECT || id == GTK_RESPONSE_ACCEPT) {
	*ret = id;
    } else if (id == GTK_RESPONSE_DELETE_EVENT) {
	*ret = GTK_RESPONSE_REJECT;
    }

    if (*ret == GTK_RESPONSE_ACCEPT) {
	struct mail_dialog *md = g_object_get_data(G_OBJECT(dlg), "md");

	err = finalize_mail_settings(md);
    }

    if (!err) {
	gtk_widget_destroy(GTK_WIDGET(dlg));
    }
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

static GtkWidget *mailer_help_button (GtkWidget *hbox,
				      void (*help_func))
{
    GtkWidget *button;

    button = gtk_button_new_from_stock(GTK_STOCK_HELP);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       button, TRUE);
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION >= 2
    gtk_button_box_set_child_non_homogeneous(GTK_BUTTON_BOX(hbox),
					     button, TRUE);
#endif
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(help_func), NULL);

    return button;
}

/* Below: notebook-style dialog with two pages, one pertaining to
   the current message and one for general "mail setup" (server
   and credentials).
*/

static int mail_to_dialog (const char *fname,
			   struct mail_info *minfo,
			   struct msg_info *msg,
			   GtkWindow *parent,
			   void (*help_func))
{
    struct mail_dialog md;
    const gchar *lbls[] = {
	N_("To:"),
	N_("From:"),
	N_("Subject:"),
	N_("Note:")
    };
    GtkWidget *tbl, *lbl;
    GtkWidget *hbox, *vbox;
    gchar *shortname = NULL;
    int i, datafile, nrows;
    gint ret = 0;

    md.dlg =
	gtk_dialog_new_with_buttons(_("gretl: send mail"),
				    parent,
				    GTK_DIALOG_DESTROY_WITH_PARENT,
				    GTK_STOCK_CANCEL,
				    GTK_RESPONSE_REJECT,
				    GTK_STOCK_OK,
				    GTK_RESPONSE_ACCEPT,
				    NULL);
    g_object_set_data(G_OBJECT(md.dlg), "md", &md);

    md.minfo = minfo;
    md.msg = msg;

    read_mail_info(minfo);
    minfo->sig = get_signature();
    minfo->want_sig = minfo->sig != NULL;

    g_signal_connect(G_OBJECT(md.dlg), "destroy",
		     G_CALLBACK(mail_dialog_quit), NULL);

    set_dialog_border_widths(md.dlg);
#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG(md.dlg), FALSE);
#endif
    gtk_window_set_position(GTK_WINDOW(md.dlg), GTK_WIN_POS_MOUSE);

    /* set up the notebook */
    md.nb = gtk_notebook_new();
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(md.dlg));
    gtk_container_add(GTK_CONTAINER(vbox), md.nb);

    /* the first page */
    hbox = gtk_hbox_new(FALSE, 5);
    border_width(hbox, 5);
    vbox = gtk_vbox_new(FALSE, 5);
    border_width(vbox, 5);
    gtk_container_add(GTK_CONTAINER(hbox), vbox);
    lbl = gtk_label_new(_("Message"));
    gtk_notebook_append_page(GTK_NOTEBOOK(md.nb), hbox, lbl);

    nrows = (minfo->sig == NULL)? 4 : 5;

    /* table within the first page */
    tbl = gtk_table_new(nrows, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_container_add(GTK_CONTAINER(vbox), tbl);

    /* what sort of file are we sending? */
    shortname = g_path_get_basename(fname);
    datafile = is_data_file(shortname);

    /* loop across data enty widgets on page 1 */
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
				       shortname);
	    } else {
		note = g_strdup_printf(_("Please find the gretl script %s attached."),
				       shortname);
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
	    md.sender_entry = w;
	} else if (i == 2) {
	    md.subj_entry = w;
	} else {
	    md.note_entry = w;
	}
    }

    g_free(shortname);

    /* conditional addition of another widget on page 1 */
    if (minfo->sig != NULL) {
	GtkWidget *w;

	w = gtk_check_button_new_with_label(_("Append signature"));
	g_signal_connect(G_OBJECT(w), "toggled", G_CALLBACK(sig_callback), minfo);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), w, 0, 2, 4, 5);
	i++;
    }

    /* setup for second page */
    hbox = gtk_hbox_new(FALSE, 5);
    border_width(hbox, 5);
    vbox = gtk_vbox_new(FALSE, 5);
    border_width(vbox, 5);
    gtk_container_add(GTK_CONTAINER(hbox), vbox);
    lbl = gtk_label_new(_("Mail setup"));
    gtk_notebook_append_page(GTK_NOTEBOOK(md.nb), hbox, lbl);

    /* table for second page */
    tbl = gtk_table_new(3, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, 0, 0, 0);

    /* SMTP server entry */
    lbl = gtk_label_new(_("SMTP server:"));
    gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
    gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, 0, 1, GTK_FILL, GTK_FILL, 0, 0);
    md.server_entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), md.server_entry, 1, 2, 0, 1);
    if (minfo->server != NULL) {
	gtk_entry_set_text(GTK_ENTRY(md.server_entry), minfo->server);
    }

    /* email userID entry */
    lbl = gtk_label_new(_("Mail username:"));
    gtk_misc_set_alignment(GTK_MISC(lbl), 1, 0.5);
    gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 1, 1, 2, GTK_FILL, GTK_FILL, 0, 0);
    md.user_entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), md.user_entry, 1, 2, 1, 2);
    if (minfo->mail_user != NULL) {
	gtk_entry_set_text(GTK_ENTRY(md.user_entry), minfo->mail_user);
    }

    /* "password required" checkbox */
    md.pass_check = gtk_check_button_new_with_label
	(_("A password is required for sending mail"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(md.pass_check),
				 minfo->pass_needed);
    gtk_table_attach(GTK_TABLE(tbl), md.pass_check, 0, 2, 2, 3,
		     GTK_FILL, GTK_FILL, 0, 0);

    /* Help button */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(md.dlg));
    mailer_help_button(hbox, help_func);

    gtk_widget_set_size_request(md.dlg, 460, -1);

    g_signal_connect(G_OBJECT(md.dlg), "response",
		     G_CALLBACK(mail_dialog_callback), &ret);
    gtk_widget_show_all(md.dlg);

    if (minfo->server == NULL || minfo->mail_user == NULL) {
	/* Mail setup not completed yet */
	gtk_notebook_set_current_page(GTK_NOTEBOOK(md.nb), 1);
    }

    gtk_main(); /* blocking, but not modal */

    return (ret == GTK_RESPONSE_ACCEPT)? 1 : 0;
}

static void toggle_store (GtkToggleButton *b, struct mail_info *minfo)
{
    gboolean store = gtk_toggle_button_get_active(b);

    if (!store) {
	minfo->store = 0;
    } else {
	GtkWidget *d;
	gint resp;

	d = gtk_message_dialog_new(NULL,
				   GTK_DIALOG_MODAL,
				   GTK_MESSAGE_QUESTION,
				   GTK_BUTTONS_YES_NO,
				   "%s",
				   _("Really store password?"));
	resp = gtk_dialog_run(GTK_DIALOG(d));
	gtk_widget_destroy(d);
	if (resp == GTK_RESPONSE_YES) {
	    minfo->store = 1;
	}
    }
}

static void cancel_password (GtkWidget *w, struct passwd_dialog *pwd)
{
    gtk_widget_destroy(pwd->dlg);
}

static void set_password (GtkWidget *w, struct passwd_dialog *pwd)
{
    const gchar *word = gtk_entry_get_text(GTK_ENTRY(pwd->entry));

    if (pwd->minfo->mail_pass == NULL ||
	strcmp(word, pwd->minfo->mail_pass)) {
	pwd->msg->changed = 1;
	g_free(pwd->minfo->mail_pass);
	pwd->minfo->mail_pass = g_strdup(word);
    }
    pwd->msg->canceled = 0;
    gtk_widget_destroy(pwd->dlg);
}

static void icon_press_callback (GtkEntry *entry,
				 GtkEntryIconPosition pos,
				 GdkEvent *event,
				 gpointer user_data)
{
    gboolean vis = gtk_entry_get_visibility(entry);

    gtk_entry_set_visibility(entry, !vis);
}

static void password_dialog (struct msg_info *msg,
			     struct mail_info *minfo)
{
    GtkWidget *hbox, *vbox;
    GtkWidget *lbl, *button;
    GdkPixbuf *pbuf = NULL;
    struct passwd_dialog pwd;

    pwd.dlg = gtk_dialog_new();
    pwd.minfo = minfo;
    pwd.msg = msg;

    /* this will be controverted by clicking OK */
    msg->canceled = 1;

    pbuf = gdk_pixbuf_new_from_xpm_data((const char **) eye_xpm);

    gtk_window_set_title(GTK_WINDOW(pwd.dlg), _("gretl: mail password"));
    set_dialog_border_widths(pwd.dlg);
    gtk_window_set_position(GTK_WINDOW(pwd.dlg), GTK_WIN_POS_MOUSE);
    g_signal_connect(G_OBJECT(pwd.dlg), "destroy",
		     G_CALLBACK(mail_dialog_quit), NULL);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(pwd.dlg));
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, 0, 0, 5);

    lbl = gtk_label_new(_("Email password:"));
    gtk_box_pack_start(GTK_BOX(hbox), lbl, 0, 0, 5);
    pwd.entry = gtk_entry_new();
    gtk_entry_set_visibility(GTK_ENTRY(pwd.entry), FALSE);
    if (minfo->mail_pass != NULL) {
	gtk_entry_set_text(GTK_ENTRY(pwd.entry), minfo->mail_pass);
    }
    gtk_entry_set_activates_default(GTK_ENTRY(pwd.entry), TRUE);
    gtk_entry_set_icon_from_pixbuf(GTK_ENTRY(pwd.entry),
				   GTK_ENTRY_ICON_SECONDARY,
				   pbuf);
    gtk_entry_set_icon_activatable(GTK_ENTRY(pwd.entry),
				   GTK_ENTRY_ICON_SECONDARY,
				   TRUE);
    g_signal_connect(G_OBJECT(pwd.entry), "icon-press",
		     G_CALLBACK(icon_press_callback), NULL);
    gtk_container_add(GTK_CONTAINER(hbox), pwd.entry);

    /* check button for storing password */
    button = gtk_check_button_new_with_label(_("Store password"));
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), button, 0, 0, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, 0, 0, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
				 minfo->store);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(toggle_store), minfo);

    /* button zone */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(pwd.dlg));

    /* Cancel button */
    button = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(cancel_password), &pwd);

    /* "OK" button */
    button = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_password), &pwd);
    gtk_widget_grab_default(button);

    gtk_widget_set_size_request(pwd.dlg, 360, -1);
    gtk_widget_show_all(pwd.dlg);

    gtk_window_set_modal(GTK_WINDOW(pwd.dlg), TRUE);
    gtk_main();
}

static void mail_infobox (const char *msg)
{
    GtkWidget *dialog;

    dialog = gtk_message_dialog_new(NULL,
				    GTK_DIALOG_DESTROY_WITH_PARENT,
				    GTK_MESSAGE_INFO,
				    GTK_BUTTONS_CLOSE,
				    "%s",
				    msg);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}

static void get_password (struct msg_info *msg,
			  struct mail_info *minfo)
{
    int old_store = minfo->store;

    password_dialog(msg, minfo);

    if (!msg->canceled) {
	if (minfo->store) {
	    if (!old_store || msg->changed) {
		write_mail_info(minfo);
	    }
	} else if (old_store) {
	    write_mail_info(minfo);
	}
    }
}

static int pack_and_mail (const char *fname,
			  struct msg_info *msg,
			  struct mail_info *minfo)
{
    char tmpfname[FILENAME_MAX];
    FILE *fpin, *fpout;
    int err = 0;

    fpin = gretl_fopen(fname, "rb");
    if (fpin == NULL) {
	perror(fname);
	err = 1;
    }

    sprintf(tmpfname, "%smpack.XXXXXX", gretl_dotdir());
    fpout = gretl_mktemp(tmpfname, "wb");
    if (fpout == NULL) {
	err = 1;
    }

    if (!err) {
	const char *ctype = is_data_file(fname) ?
	    "application/x-gretldata" : "application/x-gretlscript";

	err = mpack_encode(fpin, fname, msg->note, msg->subj,
			   msg->recip, msg->sender, ctype, fpout);
    }

    if (fpin != NULL) {
	fclose(fpin);
    }

    if (fpout != NULL) {
	fclose(fpout);
    }

    if (!err && minfo->pass_needed) {
	get_password(msg, minfo);
    }

    if (!err && !msg->canceled) {
#if 0
	fprintf(stderr, "calling curl_send_mail:\n");
	fprintf(stderr, "  sender = '%s'\n", msg->sender);
	fprintf(stderr, "  recipient = '%s'\n", msg->recip);
	fprintf(stderr, "  server = '%s'\n", minfo->server);
	fprintf(stderr, "  username = '%s'\n", minfo->mail_user);
	fprintf(stderr, "  password = (%d bytes)\n",
		(int) strlen(minfo->mail_pass));
	fprintf(stderr, "  filename = '%s'\n", tmpfname);
#endif
	err = curl_send_mail(msg->sender, msg->recip,
			     minfo->server, minfo->mail_user,
			     minfo->mail_pass, tmpfname);
    }

    gretl_remove(tmpfname);

    return err;
}

int email_file (const char *fname, GtkWindow *parent, void (*help_func))
{
    struct mail_info *minfo = NULL;
    struct msg_info *msg = NULL;
    int doit = 0;
    int err = 0;

    minfo = mail_info_new();
    msg = mail_msg_new();

    if (minfo == NULL || msg == NULL) {
	mail_msg_free(msg);
	mail_info_free(minfo);
	return E_ALLOC;
    }

    doit = mail_to_dialog(fname, minfo, msg, parent, help_func);

    if (doit) {
	err = pack_and_mail(fname, msg, minfo);
	if (!err && !msg->canceled) {
	    mail_infobox(_("Mail sent"));
	}
    }

    mail_msg_free(msg);
    mail_info_free(minfo);

    return err;
}
