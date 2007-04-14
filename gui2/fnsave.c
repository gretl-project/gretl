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
#include "datafiles.h"
#include "textbuf.h"
#include "fileselect.h"
#include "gretl_www.h"
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
    GtkWidget *ifsel;
    GtkWidget *codesel;
    GtkWidget *check;
    GtkWidget *fcheck;
    fnpkg *pkg;
    char *fname;
    char *author;
    char *version;
    char *date;
    char *pkgdesc;
    char *help;
    int pub;
    int *privlist;
    int iface;
    FuncDataReq dreq;
    float minver;
    int upload;
    int usever;
    int saveas;
};

struct login_info_ {
    GtkWidget *dlg;
    GtkWidget *login_entry;
    GtkWidget *pass_entry;
    char *login;
    char *pass;
    int canceled;
};

function_info *finfo_new (void)
{
    function_info *finfo;

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) {
	return NULL;
    }

    finfo->pkg = NULL;
    finfo->fname = NULL;
    finfo->author = NULL;
    finfo->version = NULL;
    finfo->date = NULL;
    finfo->pkgdesc = NULL;
    finfo->upload = 0;
    finfo->usever = 0;
    finfo->saveas = 0;
    finfo->iface = -1;

    finfo->help = NULL;
    finfo->pub = -1;
    finfo->privlist = NULL;
    finfo->dreq = 0;
    finfo->minver = 1.6;

    return finfo;
}

static void finfo_free (function_info *finfo)
{
    free(finfo->fname);
    free(finfo->author);
    free(finfo->version);
    free(finfo->date);
    free(finfo->pkgdesc);
    free(finfo->privlist);
    free(finfo->help);

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

static void login_finalize (GtkWidget *w, login_info *linfo)
{
    linfo->login = entry_box_get_trimmed_text(linfo->login_entry);
    linfo->pass = entry_box_get_trimmed_text(linfo->pass_entry);

    gtk_widget_destroy(linfo->dlg);
}

/* for use in File, Save dialog */

void get_default_package_name (char *fname, gpointer p)
{
    function_info *finfo = (function_info *) p;
    const char *pubname;

    *fname = '\0';
    pubname = user_function_name_by_index(finfo->pub);  

    if (pubname != NULL) {
	strcpy(fname, pubname);
	if (finfo->usever) {
	    strcat(fname, "-");
	    strcat(fname, finfo->version);
	}
	strcat(fname, ".gfn");	
    }
}

static int check_version_string (const char *s)
{
    int dc = 0;
    int err = 0;

    if (!isdigit(*s) || (*s && !isdigit(s[strlen(s) - 1]))) {
	err = 1;
    }

    while (*s && !err) {
	if (*s == '.' && ++dc > 2) {
	    err = 1;
	} else if (!isdigit(*s) && strspn(s, ".") != 1) {
	    err = 1;
	}
	s++;
    }

    return err;
}

static int maybe_revise_package_name (function_info *finfo)
{
    gchar *base = g_strdup(finfo->fname);
    const char *pubname;
    char *p;

    p = strrchr(base, SLASH);
    if (p != NULL) {
	*(p + 1) = '\0';
    }

    free(finfo->fname);
    pubname = user_function_name_by_index(finfo->pub);

    if (finfo->usever) {
	/* factor in version string */
	finfo->fname = g_strdup_printf("%s%s-%s.gfn", base, 
				       pubname, finfo->version);
    } else {
	/* no version string */
	finfo->fname = g_strdup_printf("%s%s.gfn", base, 
				       pubname);
    }
    
    g_free(base);

    return 0;
}

static void real_finfo_save (function_info *finfo)
{
    char **fields[] = {
	&finfo->author,
	&finfo->version,
	&finfo->date,
	&finfo->pkgdesc
    };
    int i, hidx = 0;
    int err = 0;

    for (i=0; i<NENTRIES && !err; i++) {
	free(*fields[i]);
	*fields[i] = entry_box_get_trimmed_text(finfo->entries[i]);
	if (*fields[i] == NULL) {
	    err = 1;
	}
    }

    if (err) {
	errbox(_("Some required information is missing"));
	return;
    }

    if (check_version_string(finfo->version)) {
	errbox(_("Invalid version string: use numbers and '.' only"));
	return;
    }

    finfo->upload = 
	gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(finfo->check));

    finfo->usever = 
	gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(finfo->fcheck));

    if (hidx >= 0) {
	char *tmp = textview_get_text(finfo->text);

	finfo->help = trim_text(tmp);
	free(tmp);
    }

    if (finfo->saveas) {
	file_selector(_("Save function package"), SAVE_FUNCTIONS, 
		      FSEL_DATA_MISC, finfo);
    } else {
	maybe_revise_package_name(finfo);
	err = save_user_functions(finfo->fname, finfo);
	if (!err) {
	    infobox(_("Saved package as %s"), finfo->fname);
	}
    }
}

static void finfo_save (GtkWidget *w, function_info *finfo)
{
    finfo->saveas = (finfo->fname == NULL);
    real_finfo_save(finfo);
}

static void finfo_destroy (GtkWidget *w, function_info *finfo)
{
    finfo_free(finfo);
}

static void login_cancel (GtkWidget *w, login_info *linfo)
{
    linfo->canceled = 1;
    gtk_widget_destroy(linfo->dlg);
}

enum {
    HIDX_INIT,
    HIDX_SWITCH
};

const char *print_today (void)
{
    static char timestr[16];
    struct tm *local;
    time_t t;

    t = time(NULL);
    local = localtime(&t);
    strftime(timestr, 15, "%Y-%m-%d", local);

    return timestr;
}

static gboolean update_iface (GtkOptionMenu *menu, 
			      function_info *finfo)
{
    int i = gtk_option_menu_get_history(menu);

    if (i == 0) {
	finfo->iface = finfo->pub;
    } else {
	finfo->iface = finfo->privlist[i];
    }

    return FALSE;
}

static void edit_code_callback (GtkWidget *w, function_info *finfo)
{
    windata_t *vwin;
    GtkWidget *orig;
    const char *funname;
    GQuark q;
    PRN *prn = NULL;

    if (finfo->iface < 0) {
	return;
    }

    funname = user_function_name_by_index(finfo->iface);
    q = g_quark_from_string(funname);

    orig = match_window_by_data(GINT_TO_POINTER(q));
    if (orig != NULL) {
	gtk_window_present(GTK_WINDOW(orig));
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    gretl_function_print_code(finfo->iface, prn);

    vwin = view_buffer(prn, 78, 350, funname,
		       EDIT_FUNC_CODE, GINT_TO_POINTER(q));

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
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 0);

    label = gtk_label_new(txt);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_widget_show(label);

    return hbox;
}

enum {
    REGULAR_BUTTON,
    CHECK_BUTTON
};

static GtkWidget *button_in_hbox (GtkWidget *w, int btype, const char *txt)
{
    GtkWidget *hbox, *button;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 0);
    if (btype == CHECK_BUTTON) {
	button = gtk_check_button_new_with_label(txt);
    } else {
	button = gtk_button_new_with_label(txt);
    }
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_widget_show(button);
    gtk_widget_show(hbox);

    return button;
}

enum {
    IFACE_PUBLIC,
    IFACE_ALL
};

static GtkWidget *interface_selector (function_info *finfo, int iface)
{
    GtkWidget *ifmenu, *menu, *tmp;
    const char *fnname;
    int i;

    finfo->iface = finfo->pub;

    ifmenu = gtk_option_menu_new();
    menu = gtk_menu_new();

    fnname = user_function_name_by_index(finfo->pub);
    tmp = gtk_menu_item_new_with_label(fnname);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), tmp);

    if (iface == IFACE_ALL && finfo->privlist != NULL) {
	for (i=1; i<=finfo->privlist[0]; i++) {
	    fnname = user_function_name_by_index(finfo->privlist[i]);
	    tmp = gtk_menu_item_new_with_label(fnname);
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), tmp);
	}
    }	

    gtk_option_menu_set_menu(GTK_OPTION_MENU(ifmenu), menu);
    gtk_widget_show_all(ifmenu);

    return ifmenu;
}

static void dreq_select (GtkOptionMenu *menu, function_info *finfo)
{
    finfo->dreq = gtk_option_menu_get_history(menu);
}

static void add_data_requirement_menu (GtkWidget *tbl, int i, 
				       function_info *finfo)
{
    const char *datareq[] = {
	N_("No special requirement"),
	N_("Time-series data"),
	N_("Quarterly or monthly data"),
	N_("Panel data")
    };
    GtkWidget *menu, *datamenu, *tmp;
    int j;

    tmp = gtk_label_new(_("Data requirement"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
    gtk_widget_show(tmp);

    datamenu = gtk_option_menu_new();
    menu = gtk_menu_new();
    for (j=0; j<=FN_NEEDS_PANEL; j++) {
	tmp = gtk_menu_item_new_with_label(_(datareq[j]));
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), tmp);
    }
    gtk_option_menu_set_menu(GTK_OPTION_MENU(datamenu), menu);
    gtk_option_menu_set_history(GTK_OPTION_MENU(datamenu), finfo->dreq);

    tmp = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tmp), datamenu, FALSE, FALSE, 0);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 1, 2, i, i+1);
    gtk_widget_show_all(tmp);

    g_signal_connect(G_OBJECT(datamenu), "changed",
		     G_CALLBACK(dreq_select), finfo);
}

static void get_maj_min_pl (float minver, int *maj, int *min, int *pl)
{
    char vstr[5], minstr[2], plstr[2];

    gretl_push_c_numeric_locale();
    sprintf(vstr, "%.2f", (double) minver);
    gretl_pop_c_numeric_locale();

    sscanf(vstr, "%d.%1s%1s", maj, minstr, plstr);
    *min = atoi(minstr);
    *pl = atoi(plstr);
}

static void adjust_minver (GtkWidget *w, function_info *finfo)
{
    int val = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int lev = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "level"));
    int maj, min, pl;

    get_maj_min_pl(finfo->minver, &maj, &min, &pl);

    if (lev == 1) {
	finfo->minver = (float) val + min / 10.0 + pl / 100.0;
    } else if (lev == 2) {
	finfo->minver = (float) maj + val / 10.0 + pl / 100.0;
    } else if (lev == 3) {
	finfo->minver = (float) maj + min / 10.0 + val / 100.0;
    }
}

static void add_minver_selector (GtkWidget *tbl, int i, 
				 function_info *finfo)
{
    GtkWidget *tmp, *spin, *hbox;
    int maj, min, pl;

    get_maj_min_pl(finfo->minver, &maj, &min, &pl);

    tmp = gtk_label_new(_("Minimum gretl version"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
    gtk_widget_show(tmp);

    hbox = gtk_hbox_new(FALSE, 0);

    spin = gtk_spin_button_new_with_range(1, 3, 1);
    if (maj > 1) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), (double) maj);
    }
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 2);
    g_object_set_data(G_OBJECT(spin), "level", GINT_TO_POINTER(1));
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(adjust_minver), finfo);
    tmp = gtk_label_new(".");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 2);

    spin = gtk_spin_button_new_with_range(0, 9, 1);
    if (min > 0) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), (double) min);
    }
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 2);
    g_object_set_data(G_OBJECT(spin), "level", GINT_TO_POINTER(2));
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(adjust_minver), finfo);
    tmp = gtk_label_new(".");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 2);

    spin = gtk_spin_button_new_with_range(0, 9, 1);
    if (pl > 0) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), (double) pl);
    }
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 2);
    g_object_set_data(G_OBJECT(spin), "level", GINT_TO_POINTER(3));
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(adjust_minver), finfo);

    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 1, 2, i, i+1);
    gtk_widget_show_all(hbox);
}

static GtkWidget *editable_text_box (void)
{
    GtkTextBuffer *tbuf = gretl_text_buf_new();
    GtkWidget *w = gtk_text_view_new_with_buffer(tbuf);

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(w), 4);
    gtk_widget_modify_font(GTK_WIDGET(w), fixed_font);
    gtk_text_view_set_editable(GTK_TEXT_VIEW(w), TRUE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w), TRUE);

    return w;
}

static const gchar *get_user_string (void)
{
    const gchar *name;

    name = g_get_real_name();
    if (name == NULL) {
	name = g_get_user_name();
    }

    return name;
}

static void finfo_dialog (function_info *finfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *vbox, *hbox;
    const char *entry_labels[] = {
	N_("Author"),
	N_("Version"),
	N_("Date (YYYY-MM-DD)"),
	N_("Package description")
    };
    char *entry_texts[] = {
	finfo->author,
	finfo->version,
	finfo->date,
	finfo->pkgdesc
    };
    const char *fnname;
    int focused = 0;
    const char *hlp;
    gchar *ltxt;
    int i;

    finfo->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    g_signal_connect(G_OBJECT(finfo->dlg), "destroy", 
		     G_CALLBACK(finfo_destroy), finfo);

    gtk_window_set_title(GTK_WINDOW(finfo->dlg), 
			 _("gretl: function package editor")); 
    gtk_window_set_default_size(GTK_WINDOW(finfo->dlg), 640, 500);

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);
    gtk_container_add(GTK_CONTAINER(finfo->dlg), vbox);
    gtk_widget_show(vbox);
			 
    tbl = gtk_table_new(NENTRIES + 1, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 4);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 5);

    for (i=0; i<NENTRIES; i++) {
	GtkWidget *entry;

	label = gtk_label_new(_(entry_labels[i]));
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, i, i+1);
	gtk_widget_show(label);

	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 40);
	gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	gtk_widget_show(entry); 

	finfo->entries[i] = entry;

	if (entry_texts[i] != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(entry), entry_texts[i]);
	} else if (i == 0) {
	    const gchar *s = get_user_string();

	    if (s != NULL) {
		gtk_entry_set_text(GTK_ENTRY(entry), s);
	    }
	} else if (i == 1) {
	    gtk_entry_set_text(GTK_ENTRY(entry), "1.0");
	} else if (i == 2) {
	    gtk_entry_set_text(GTK_ENTRY(entry), print_today());
	}

	if (i == 0 && entry_texts[i] == NULL) {
	    gtk_widget_grab_focus(entry);
	    focused = 1;
	} else if (i == 1 && !focused) {
	    gtk_widget_grab_focus(entry);
	}
    }

    add_minver_selector(tbl, i++, finfo);
    add_data_requirement_menu(tbl, i, finfo);
    gtk_widget_show(tbl);

    fnname = user_function_name_by_index(finfo->pub);
    ltxt = g_strdup_printf(_("Help text for %s:"), fnname);
    hbox = label_hbox(vbox, ltxt);
    gtk_widget_show(hbox);
    g_free(ltxt);

    finfo->text = editable_text_box();
    text_table_setup(vbox, finfo->text);

    gretl_function_get_info(finfo->pub, "help", &hlp);
    textview_set_text(finfo->text, hlp);

    /* edit code button, possibly with selector */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    button = gtk_button_new_with_label(_("Edit function code"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_code_callback), finfo);

    if (finfo->privlist != NULL) {
	finfo->codesel = interface_selector(finfo, IFACE_ALL);
	gtk_box_pack_start(GTK_BOX(hbox), finfo->codesel, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(GTK_OPTION_MENU(finfo->codesel)), "changed",
			 G_CALLBACK(update_iface), finfo);
    } else {
	finfo->iface = finfo->pub;
    }

    gtk_widget_show_all(hbox);

    /* check box for upload option */
    finfo->check = button_in_hbox(vbox, CHECK_BUTTON, 
				  _("Upload package to server on save"));

    /* check box for use version number in filename */
    finfo->fcheck = button_in_hbox(vbox, CHECK_BUTTON, 
				   _("Include version in file name"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(finfo->fcheck),
				 finfo->usever);    

    /* control button area */
    hbox = gtk_hbutton_box_new();
    gtk_button_box_set_layout(GTK_BUTTON_BOX(hbox), GTK_BUTTONBOX_END);
    gtk_button_box_set_spacing(GTK_BUTTON_BOX(hbox), 10);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* Save button */
    button = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(finfo_save), finfo);
    gtk_widget_grab_default(button);

    /* Close button */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(delete_widget), finfo->dlg);

    gtk_widget_show_all(hbox);

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
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(linfo->dlg)->vbox), tbl, 
		       FALSE, FALSE, 5);

    for (i=0; i<2; i++) {
	char *src = (i == 0)? linfo->login : linfo->pass;
	GtkWidget *entry;

	label = gtk_label_new((i == 0)? _("Login") : _("Password"));
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, i, i+1,
			 GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL,
			 5, 5);
	gtk_widget_show(label);

	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 34);
	gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
	gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
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

    /* control button area */

    hbox = GTK_DIALOG(linfo->dlg)->action_area;

    /* Cancel */
    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(login_cancel), linfo);
    gtk_widget_show(button);

    /* OK */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(login_finalize), linfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* Website */
    button = gtk_button_new_with_label("Website");
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       button, TRUE);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(web_get_login), NULL);
    gtk_widget_show(button);

    gtk_widget_show(linfo->dlg);
}

#define url_reserved(c) (strchr(";/?:@&=+$,<>%#\t\r\n\v\0", c) != NULL)

static int count_specials (const char *s)
{
    int n = 0;

    while (*s) {
	if (url_reserved(*s) || !isprint(*s)) {
	    n++;
	}
	s++;
    }

    return n;
}

static char *url_encode_string (const char *s)
{
    char *encstr, *p;
    int n;

    if (s == NULL) {
	return NULL;
    }

    n = count_specials(s);
    if (n == 0) {
	return gretl_strdup(s);
    }

    encstr = malloc(strlen(s) + n * 2 + 1);

    if (encstr != NULL) {
	p = encstr;
	while (*s) {
	    if (*s == ' ') {
		*p++ = '+';
	    } else if (url_reserved(*s) || !isprint(*s)) {
		sprintf(p, "%%%.2X", *s);
		p += 3;
	    } else {
		*p++ = *s;
	    } 
	    s++;
	}
	*p = '\0';
    }

    return encstr;
}

static void do_upload (const char *fname)
{
    char *ulogin = NULL;
    char *upass = NULL;
    char *ufname = NULL;
    char *ubuf = NULL;
    char *buf = NULL;
    char *retbuf = NULL;
    login_info linfo;
    GdkDisplay *disp;
    GdkCursor *cursor;
    GdkWindow *w1;
    gint x, y;
    int err = 0;

    login_dialog(&linfo);

    if (linfo.canceled) {
	linfo_free(&linfo);
	return;
    }

    /* set waiting cursor */
    disp = gdk_display_get_default();
    cursor = gdk_cursor_new(GDK_WATCH);
    w1 = gdk_display_get_window_at_pointer(disp, &x, &y);
    gdk_window_set_cursor(w1, cursor);
    gdk_display_sync(disp);
    gdk_cursor_unref(cursor);

    g_file_get_contents(fname, &buf, NULL, NULL);

    if (buf == NULL) {
	err = E_ALLOC;
    } else {  
	ubuf = url_encode_string(buf);
	ulogin = url_encode_string(linfo.login);
	upass = url_encode_string(linfo.pass);
	ufname = url_encode_string(path_last_element(fname));
	
	if (ubuf == NULL || ulogin == NULL || upass == NULL || ufname == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = upload_function_package(ulogin, upass, ufname, ubuf,
				      &retbuf);
    }

    /* reset default cursor */
    gdk_window_set_cursor(w1, NULL);

    if (err) {
	gui_errmsg(err);
    } else if (retbuf != NULL && *retbuf != '\0') {
	infobox(retbuf);
    }

    free(ulogin);
    free(upass);
    free(ufname);
    free(ubuf);
    free(buf);
    free(retbuf);

    linfo_free(&linfo);
}

static int 
fnpkg_check_filename (function_info *finfo, const char *fname)
{
    const char *p = strrchr(fname, SLASH);
    char funname[FN_NAMELEN];
    int i, n, err = 0;

    if (p == NULL) {
	p = fname;
    } else {
	p++;
    }

    n = 0;
    for (i=0; p[i] != '\0'; i++) {
	if (p[i] == '-' || !strcmp(p + i, ".gfn")) {
	    break;
	}
	n++;
    }

    if (n >= FN_NAMELEN) {
	err = 1;
    } else {
	const char *iname = user_function_name_by_index(finfo->pub);

	*funname = '\0';
	strncat(funname, p, n);
	if (strcmp(funname, iname)) {
	    err = 1;
	}
    }

    if (err) {
	errbox(_("The package filename must match the name of\n"
		 "the package's public interface"));
    }

    return err;
}

int save_user_functions (const char *fname, gpointer p)
{
    function_info *finfo = p;
    int i, err;

    /* sync/check filename with functions editor */
    if (finfo->fname == NULL) {
	err = fnpkg_check_filename(finfo, fname);
	if (err) {
	    return 1;
	}
	finfo->fname = g_strdup(fname);
    } 

    if (finfo->privlist != NULL) {
	for (i=1; i<=finfo->privlist[0]; i++) {
	    gretl_function_set_private(finfo->privlist[i], TRUE);
	}
    }

    gretl_function_set_info(finfo->pub, finfo->help);
    gretl_function_set_private(finfo->pub, FALSE);

#if 0
    fprintf(stderr, "author='%s'\n", finfo->author);
    fprintf(stderr, "version='%s'\n", finfo->version);
    fprintf(stderr, "date='%s'\n", finfo->date);
    fprintf(stderr, "pkgdesc='%s'\n", finfo->pkgdesc);
    fprintf(stderr, "finfo->pub = %d\n", finfo->pub);
    printlist(finfo->privlist, "finfo->privlist");
    fprintf(stderr, "dreq=%d\n", finfo->dreq);
    fprintf(stderr, "minver=%.2f\n", (double) finfo->minver);
#endif
		
    err = write_function_package(finfo->pkg,
				 fname,
				 finfo->pub, 
				 finfo->privlist,
				 finfo->author,
				 finfo->version,
				 finfo->date,
				 finfo->pkgdesc,
				 finfo->dreq,
				 finfo->minver);

    if (err) {
	gui_errmsg(err);
    } else {
	maybe_update_func_files_window(1);
	if (finfo->upload) {
	    do_upload(fname);
	}
    }

    return err;
}

/* called from function selection dialog: a set of functions has been
   selected and now we need to add info on author, version, etc.
*/

void prepare_functions_save (void)
{
    function_info *finfo;
    int *list = NULL;
    int i;

    if (storelist == NULL) {
	return;
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    list = gretl_list_from_string(storelist);
    if (list == NULL) {
	nomem();
	free(finfo);
	return;
    }

    finfo->pub = list[1];

    if (list[0] > 1) {
	finfo->privlist = gretl_list_new(list[0] - 1);
	if (finfo->privlist == NULL) {
	    nomem();
	    free(finfo);
	    free(list);
	    return;
	} else {
	    for (i=1; i<=finfo->privlist[0]; i++) {
		finfo->privlist[i] = list[i+1];
	    }
	    free(list);
	}
    } else {
	finfo->privlist = NULL;
    }

    /* Call dialog to do the actual editing */
    finfo_dialog(finfo);
}

void edit_function_package (const char *fname, int *loaderr)
{
    function_info *finfo;
    const char *p;
    int err = 0;

    if (!function_package_is_loaded(fname)) {
	err = load_user_function_file(fname);
	if (err) {
	    fprintf(stderr, "load_user_function_file: failed on %s\n", fname);
	    errbox(_("Couldn't open %s"), fname);
	    *loaderr = 1;
	    return;
	}
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    err = function_package_get_info(fname,
				    &finfo->pkg,
				    &finfo->pub,
				    &finfo->privlist,
				    &finfo->author,
				    &finfo->version,
				    &finfo->date,
				    &finfo->pkgdesc,
				    &finfo->dreq,
				    &finfo->minver);

    if (err) {
	fprintf(stderr, "function_package_get_info: failed on %s\n", fname);
	errbox("Couldn't get function package information");
	finfo_free(finfo);
	return;
    }

    p = strrchr(fname, SLASH);
    if (p == NULL) {
	p = fname;
    } else {
	p++;
    }

    if (strchr(p, '-')) {
	finfo->usever = 1;
    }

    finfo->fname = g_strdup(fname);

    finfo_dialog(finfo);
}
