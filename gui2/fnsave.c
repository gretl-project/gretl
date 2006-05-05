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

#include "gretl.h"
#include "dlgutils.h"
#include "textbuf.h"
#include "fileselect.h"
#include "webget.h"
#include "fnsave.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#include "gretl_func.h"

#define NENTRIES 4

typedef struct function_info_ function_info;
typedef struct login_info_ login_info;

struct function_info_ {
    GtkWidget *dlg;
    GtkWidget *entries[NENTRIES];
    GtkWidget *text;
    GtkWidget *combo;
    GtkWidget *codesel;
    GtkWidget *check;
    char *author;
    char *version;
    char *date;
    char *pkgdesc;
    char **help;
    int *publist;
    int *privlist;
    int n_public;
    int iface;
    int upload;
    int canceled;
};

struct login_info_ {
    GtkWidget *dlg;
    GtkWidget *login_entry;
    GtkWidget *pass_entry;
    char *login;
    char *pass;
    int canceled;
};

static const char *fnsave_filename;

static void set_fnsave_filename (const char *fname)
{
    fnsave_filename = fname;
}

const char *get_fnsave_filename (void)
{
    return fnsave_filename;
}

function_info *finfo_new (void)
{
    function_info *finfo;

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) {
	return NULL;
    }

    finfo->author = NULL;
    finfo->version = NULL;
    finfo->date = NULL;
    finfo->pkgdesc = NULL;
    finfo->upload = 0;
    finfo->canceled = 0;
    finfo->iface = -1;

    finfo->n_public = 0;
    finfo->help = NULL;
    finfo->publist = NULL;
    finfo->privlist = NULL;

    return finfo;
}

static int finfo_init (function_info *finfo)
{
    finfo->n_public = finfo->publist[0];

    finfo->help = create_strings_array(finfo->n_public);
    if (finfo->help == NULL) {
	errbox(_("Out of memory!"));
	finfo->canceled = 1;
	return E_ALLOC;
    }

    return 0;
}

static void finfo_free (function_info *finfo)
{
    free(finfo->author);
    free(finfo->version);
    free(finfo->date);
    free(finfo->pkgdesc);
    free(finfo->publist);
    free(finfo->privlist);

    free_strings_array(finfo->help, finfo->n_public);

    free(finfo);
}

static void login_init_or_free (login_info *linfo, int freeit)
{
    static char *login;
    static char *pass;

    if (freeit) {
	if (!linfo->canceled) {
	    free(login);
	    free(pass);
	    login = g_strdup(linfo->login);
	    pass = g_strdup(linfo->pass);
	}
	free(linfo->login);
	free(linfo->pass);
    } else {
	linfo->login = (login == NULL)? NULL : g_strdup(login);
	linfo->pass = (pass == NULL)? NULL : g_strdup(pass);
	linfo->canceled = 0;
    }
}

static void login_init (login_info *linfo)
{
    login_init_or_free(linfo, 0);
}

static void linfo_free (login_info *linfo)
{
    login_init_or_free(linfo, 1);
}

static char *trim_text (const char *s)
{
    char *ret = NULL;
    int i, len;

    while (isspace(*s)) s++;
    if (*s == '\0') return NULL;

    len = strlen(s);
    for (i=len-1; i>0; i--) {
	if (!isspace(s[i])) break;
	len--;
    }

    if (len > 0) {
	ret = g_strndup(s, len);
    }

    return ret;
}

static int help_text_index (function_info *finfo)
{
    const char *fname;
    int i, idx;

    fname = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(finfo->combo)->entry));
    idx = user_function_index_by_name(fname);

    for (i=0; i<finfo->n_public; i++) {
	if (idx == finfo->publist[i+1]) {
	    return i;
	}
    }

    return -1;
}

static void login_finalize (GtkWidget *w, login_info *linfo)
{
    const gchar *txt;

    txt = gtk_entry_get_text(GTK_ENTRY(linfo->login_entry));
    if (txt != NULL && *txt != '\0') {
	linfo->login = trim_text(txt);
    }

    txt = gtk_entry_get_text(GTK_ENTRY(linfo->pass_entry));
    if (txt != NULL && *txt != '\0') {
	linfo->pass = trim_text(txt);
    }

    gtk_widget_destroy(linfo->dlg);
}

static void finfo_finalize (GtkWidget *w, function_info *finfo)
{
    char *fields[] = {
	finfo->author,
	finfo->version,
	finfo->date,
	finfo->pkgdesc
    };
    const gchar *txt;
    int i, hidx = 0;
    int err = 0;

    for (i=0; i<NENTRIES && !err; i++) {
	free(fields[i]);
	fields[i] = NULL;
	txt = gtk_entry_get_text(GTK_ENTRY(finfo->entries[i]));
	if (txt != NULL && *txt != '\0') {
	    fields[i] = trim_text(txt);
	} else {
	    err = 1;
	}
    }

    if (err) {
	errbox(_("Some required information is missing"));
	return;
    }

    finfo->upload = 
	gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(finfo->check));

    if (finfo->n_public > 1) {
	hidx = help_text_index(finfo);
    } 

    if (hidx >= 0) {
	char *tmp = textview_get_text(finfo->text);

	finfo->help[hidx] = trim_text(tmp);
	free(tmp);
    }

    gtk_widget_destroy(finfo->dlg);
}

static void login_cancel (GtkWidget *w, login_info *linfo)
{
    linfo->canceled = 1;
    gtk_widget_destroy(linfo->dlg);
}

static void finfo_cancel (GtkWidget *w, function_info *finfo)
{
    finfo->canceled = 1;
    gtk_widget_destroy(finfo->dlg);
}

static void finfo_delete (GtkWidget *w, function_info *finfo)
{
    finfo->canceled = 1;
}

enum {
    HIDX_INIT,
    HIDX_SWITCH
};

static void 
set_dialog_info_from_fn (function_info *finfo, int idx, int code)
{
    const char *attrib = NULL;
    const char *keys[] = {
	"author",
	"version",
	"date",
	"pkgdesc"
    };
    const char *etxt;

    static int old_hidx;
    int i, new_hidx;

    if (code == HIDX_INIT) {
	old_hidx = new_hidx = 0;
    } else {
	new_hidx = help_text_index(finfo);
    }

    for (i=0; i<NENTRIES; i++) {
	etxt = gtk_entry_get_text(GTK_ENTRY(finfo->entries[i]));
	if (*etxt == '\0') {
	    gretl_function_get_info(idx, keys[i], &attrib);
	    if (attrib != NULL) {
		etxt = gtk_entry_get_text(GTK_ENTRY(finfo->entries[i]));
	    }
	}
    }

    if (new_hidx != old_hidx) {
	/* we're switching the "active" interface, so save the 
	   help text for the previous interface */
	char *old_help = textview_get_text(finfo->text);

	free(finfo->help[old_hidx]);
	finfo->help[old_hidx] = old_help;
    }

    if (code == HIDX_INIT || new_hidx != old_hidx) {
	/* initializing or switching: insert new help text */
	const char *new_help;

	gretl_function_get_info(idx, "help", &new_help);
	textview_set_text(finfo->text, new_help);
    }

    old_hidx = new_hidx;
}

static gboolean update_public (GtkEditable *entry, 
			       function_info *finfo)
{
    const char *fnname;
    int idx;

    fnname = gtk_entry_get_text(GTK_ENTRY(entry));

    if (fnname != NULL && *fnname != '\0') {
	idx = user_function_index_by_name(fnname);
	if (idx >= 0) {
	    set_dialog_info_from_fn(finfo, idx, HIDX_SWITCH);
	}
    }

    return FALSE;
}

static gboolean update_iface (GtkEditable *entry, 
			      function_info *finfo)
{
    const char *fnname;
    int idx;

    fnname = gtk_entry_get_text(GTK_ENTRY(entry));

    if (fnname != NULL && *fnname != '\0') {
	idx = user_function_index_by_name(fnname);
	if (idx >= 0) {
	    finfo->iface = idx;
	}
    }

    return FALSE;
}

static void edit_code_callback (GtkWidget *w, function_info *finfo)
{
    
    windata_t *vwin;
    PRN *prn = NULL;

    if (finfo->iface < 0) {
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    gretl_function_print_code(finfo->iface, prn);

    vwin = view_buffer(prn, 78, 350, 
		       user_function_name_by_index(finfo->iface),
		       EDIT_FUNC_CODE, NULL);

    if (vwin != NULL) {
	build_path(vwin->fname, paths.userdir, "pkgedit", NULL);
	gretl_tempname(vwin->fname);
	g_object_set_data(G_OBJECT(vwin->w), "iface", 
			  GINT_TO_POINTER(finfo->iface));
    }
}

static GtkWidget *label_hbox (GtkWidget *w, const char *txt)
{
    GtkWidget *hbox, *label;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 5);

    label = gtk_label_new(txt);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_widget_show(label);

    return hbox;
}

enum {
    REGULAR_BUTTON,
    CHECK_BUTTON
};

static GtkWidget *button_in_hbox (GtkWidget *w, int btype, const char *txt,
				  GtkWidget **phbox)
{
    GtkWidget *hbox, *button;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 5);
    if (btype == CHECK_BUTTON) {
	button = gtk_check_button_new_with_label(txt);
    } else {
	button = gtk_button_new_with_label(txt);
    }
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_widget_show(button);
    gtk_widget_show(hbox);

    if (phbox != NULL) {
	*phbox = hbox;
    }

    return button;
}

enum {
    IFACE_PUBLIC,
    IFACE_ALL
};

static GtkWidget *interface_selector (function_info *finfo, int iface)
{
    GList *fn_list = NULL;
    GtkWidget *combo;
    const char *fnname;
    int i;

    if (iface == IFACE_ALL) {
	if (finfo->privlist != NULL) {
	    finfo->iface = finfo->privlist[1];
	    for (i=1; i<=finfo->privlist[0]; i++) {
		fnname = user_function_name_by_index(finfo->privlist[i]);
		fn_list = g_list_append(fn_list, (gpointer) fnname);
	    }
	} else {
	    finfo->iface = finfo->publist[1];
	}
    }	

    for (i=1; i<=finfo->publist[0]; i++) {
	fnname = user_function_name_by_index(finfo->publist[i]);
	fn_list = g_list_append(fn_list, (gpointer) fnname);
    }

    combo = gtk_combo_new();
    gtk_combo_set_popdown_strings(GTK_COMBO(combo), fn_list); 
    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(combo)->entry), FALSE);
    gtk_widget_show(combo);
    g_list_free(fn_list);

    return combo;
}

static void finfo_dialog (function_info *finfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *hbox;
    const char *entry_labels[] = {
	N_("Author"),
	N_("Version"),
	N_("Date"),
	N_("Package description")
    };
    char *entry_texts[] = {
	finfo->author,
	finfo->version,
	finfo->date,
	finfo->pkgdesc
    };
    const char *fnname;
    int i;

    if (finfo_init(finfo)) {
	return;
    }

    finfo->dlg = gretl_dialog_new(_("gretl: function package editor"), NULL, 
				  GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    /* FIXME want label at top of dialog? */

    g_signal_connect(G_OBJECT(finfo->dlg), "delete_event",
		     G_CALLBACK(finfo_delete), finfo);

    tbl = gtk_table_new(NENTRIES, 2, FALSE);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox),
		       tbl, FALSE, FALSE, 5);

    for (i=0; i<NENTRIES; i++) {
	GtkWidget *entry;

	label = gtk_label_new(_(entry_labels[i]));
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, i, i+1);
	gtk_widget_show(label);

	entry = gtk_entry_new();
#ifndef OLD_GTK
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 40);
#endif
	gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	gtk_widget_show(entry); 

	finfo->entries[i] = entry;

	if (entry_texts[i] != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(entry), entry_texts[i]);
	}
    }

    gtk_widget_show(tbl);

    if (finfo->n_public > 1) {
	/* drop-down selector for public interfaces */
	hbox = label_hbox(GTK_DIALOG(finfo->dlg)->vbox, _("Help text for"));
	finfo->combo = interface_selector(finfo, IFACE_PUBLIC);
	gtk_box_pack_start(GTK_BOX(hbox), finfo->combo, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(GTK_COMBO(finfo->combo)->entry), "changed",
			 G_CALLBACK(update_public), finfo);
	gtk_widget_show(hbox);
    } else {
	/* only one public interface */
	gchar *ltxt;

	fnname = user_function_name_by_index(finfo->publist[1]);
	ltxt = g_strdup_printf("Help text for %s:", fnname);
	hbox = label_hbox(GTK_DIALOG(finfo->dlg)->vbox, ltxt);
	gtk_widget_show(hbox);
	g_free(ltxt);
    }

    finfo->text = create_text(finfo->dlg, -1, -1, TRUE);
#ifdef OLD_GTK
    tbl = text_table_setup(GTK_DIALOG(finfo->dlg)->vbox, finfo->text);
    gtk_widget_set_usize(finfo->dlg, 640, 440);
#else
    text_table_setup(GTK_DIALOG(finfo->dlg)->vbox, finfo->text);
    gtk_window_set_default_size(GTK_WINDOW(finfo->dlg), 640, 480);
#endif

    set_dialog_info_from_fn(finfo, finfo->publist[1], HIDX_INIT);

    /* button for editing the actual code */
    button = button_in_hbox(GTK_DIALOG(finfo->dlg)->vbox,
			    REGULAR_BUTTON, "Edit function code", &hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_code_callback), finfo);

    /* with selector if there's more than one function in package */
    if (finfo->publist[0] > 1 || finfo->privlist != NULL) {
	finfo->codesel = interface_selector(finfo, IFACE_ALL);
	gtk_box_pack_start(GTK_BOX(hbox), finfo->codesel, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(GTK_COMBO(finfo->codesel)->entry), "changed",
			 G_CALLBACK(update_iface), finfo);
    } else {
	finfo->iface = finfo->publist[1];
    }

    /* check box for upload option */
    finfo->check = button_in_hbox(GTK_DIALOG(finfo->dlg)->vbox, CHECK_BUTTON, 
				  "Upload package to server on save",
				  NULL);

    /* Create the "OK" button */
    button = ok_button(GTK_DIALOG(finfo->dlg)->action_area);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(finfo_finalize), finfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* And a Cancel button */
    button = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(finfo_cancel), finfo);
    gtk_widget_show(button);

    gtk_widget_show(finfo->dlg);
}

static void web_get_login (GtkWidget *w, gpointer p)
{
    browser_open("http://ricardo.ecn.wfu.edu/gretl/apply/");
}

static void login_dialog (login_info *linfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *hbox;
    int i;

    login_init(linfo);

    linfo->dlg = 
	gretl_dialog_new(_("gretl: upload"), NULL, GRETL_DLG_BLOCK);

    hbox = label_hbox(GTK_DIALOG(linfo->dlg)->vbox, _("Upload function package"));
    gtk_widget_show(hbox);

    tbl = gtk_table_new(2, 2, FALSE);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(linfo->dlg)->vbox), tbl, FALSE, FALSE, 5);

    for (i=0; i<2; i++) {
	char *src = (i == 0)? linfo->login : linfo->pass;
	GtkWidget *entry;

	label = gtk_label_new((i == 0)? _("Login") : _("Password"));
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, i, i+1,
			 GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL,
			 5, 5);
	gtk_widget_show(label);

	entry = gtk_entry_new();
#ifndef OLD_GTK
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 34);
#endif
	gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	if (src != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(entry), src);
	}
	gtk_widget_show(entry); 

	if (i == 0) {
	    linfo->login_entry = entry;
	} else {
	    gtk_entry_set_visibility(GTK_ENTRY(entry), FALSE);
	    linfo->pass_entry = entry;
	}
    }

    gtk_widget_show(tbl);

    hbox = label_hbox(GTK_DIALOG(linfo->dlg)->vbox, 
		      _("If you don't have a login to the gretl server\n"
			"please see http://ricardo.ecn.wfu.edu/gretl/apply/.\n"
			"The 'Website' button below should open this page\n"
			"in your web browser."));
    gtk_widget_show(hbox);


    button = ok_button(GTK_DIALOG(linfo->dlg)->action_area);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(login_finalize), linfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    button = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(linfo->dlg)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(login_cancel), linfo);
    gtk_widget_show(button);

    button = gtk_button_new_with_label("Website");
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(linfo->dlg)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(web_get_login), NULL);
    gtk_widget_show(button);

    gtk_widget_show(linfo->dlg);
}

static void do_upload (const char *fname)
{
    login_info linfo;
    char errbuf[128];

    login_dialog(&linfo);

    if (!linfo.canceled) {
	int err = upload_function_package(linfo.login,
					  linfo.pass,
					  fname,
					  errbuf);
	if (err) {
	    errbox(errbuf);
	}
    }

    linfo_free(&linfo);
}

void save_user_functions (const char *fname, gpointer p)
{
    function_info *finfo = p;
    int i, err;

    if (finfo->privlist != NULL) {
	for (i=1; i<=finfo->privlist[0]; i++) {
	    gretl_function_set_private(finfo->privlist[i], TRUE);
	}
    }

    for (i=1; i<=finfo->publist[0]; i++) {
	gretl_function_set_info(finfo->publist[i], finfo->help[i-1]);
	gretl_function_set_private(finfo->publist[i], FALSE);
    }
		
    err = write_function_package(fname,
				 finfo->privlist, 
				 finfo->publist,
				 finfo->author,
				 finfo->version,
				 finfo->date,
				 finfo->pkgdesc);

    if (err) {
	gui_errmsg(err);
    } else if (finfo->upload) {
	do_upload(fname);
    }

    finfo_free(finfo);    
}

/* called from function selection dialog: a set of functions has been
   selected and now we need to add info on author, version, etc.
*/

void prepare_functions_save (void)
{
    function_info *finfo;
    int *list = NULL;

    if (storelist == NULL) {
	return;
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    list = gretl_list_from_string(storelist);
    if (list == NULL) {
	errbox(_("Out of memory!"));
	free(finfo);
	return;
    }

    if (gretl_list_has_separator(list)) {
	if (gretl_list_split_on_separator(list, &finfo->privlist, 
					  &finfo->publist)) {
	    errbox(_("Out of memory!"));
	    free(finfo);
	    free(list);
	    return;
	} else {
	    free(list);
	}
    } else {
	finfo->publist = list;
	finfo->privlist = NULL;
    }

    finfo_dialog(finfo);

    if (finfo->canceled) {
	finfo_free(finfo);
    } else {
	set_fnsave_filename(NULL);
	file_selector(_("Save function package"), SAVE_FUNCTIONS, 
		      FSEL_DATA_MISC, finfo);
    }
}

void edit_function_package (const char *fname)
{
    function_info *finfo;
    int err = 0;

    if (!user_function_file_is_loaded(fname)) {
	err = load_user_function_file(fname);
	if (err) {
	    errbox(_("Couldn't open %s"), fname);
	    return;
	}
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    err = function_package_get_info(fname,
				    &finfo->privlist,
				    &finfo->publist,
				    &finfo->author,
				    &finfo->version,
				    &finfo->date,
				    &finfo->pkgdesc);

    if (err) {
	errbox("Couldn't get function package information");
	finfo_free(finfo);
	return;
    }

    set_fnsave_filename(fname);

    finfo_dialog(finfo);

    if (finfo->canceled) {
	finfo_free(finfo);
    } else {
	file_selector(_("Save function package"), SAVE_FUNCTIONS, 
		      FSEL_DATA_MISC, finfo);
    }    
}


