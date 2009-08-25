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
#include "dlgutils.h"
#include "datafiles.h"
#include "textbuf.h"
#include "fileselect.h"
#include "gretl_www.h"
#include "winstack.h"
#include "fnsave.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#include "gretl_func.h"

#define NENTRIES 4

typedef struct function_info_ function_info;
typedef struct login_info_ login_info;

/* package-editing dialog and associated info */

struct function_info_ {
    GtkWidget *dlg;        /* editing dialog box */
    GtkWidget *entries[NENTRIES]; /* entry boxes for author, etc. */
    GtkWidget *text;       /* text box for help */
    GtkWidget *codesel;    /* code-editing selector */
    GtkWidget *save;       /* Save button in dialog */
    GtkWidget *popup;      /* popup menu */
    windata_t *samplewin;  /* window for editing sample script */
    fnpkg *pkg;            /* pointer to package being edited */
    gchar *fname;          /* package filename */
    gchar *author;         /* package author */
    gchar *version;        /* package version number */
    gchar *date;           /* package last-revised date */
    gchar *pkgdesc;        /* package description */
    gchar *sample;         /* sample script for package */
    gchar *help;           /* package help text */
    char **pubnames;       /* names of public functions */
    char **privnames;      /* name sof private functions */
    int n_pub;             /* number of public functions */
    int n_priv;            /* number of private functions */
    char *active;          /* name of 'active' function */
    FuncDataReq dreq;      /* data requirement of package */
    int minver;            /* minimum gretl version, package */
    gboolean upload;       /* upload to server on save? */
    gboolean modified;     /* anything changed in package? */
};

/* info relating to login to server for upload */

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
    finfo->sample = NULL;

    finfo->upload = FALSE;
    finfo->modified = FALSE;

    finfo->active = NULL;
    finfo->samplewin = NULL;
    finfo->popup = NULL;

    finfo->help = NULL;
    finfo->pubnames = NULL;
    finfo->privnames = NULL;

    finfo->n_pub = 0;
    finfo->n_priv = 0;

    finfo->dreq = 0;
    finfo->minver = 10804;

    return finfo;
}

static char *funname_from_filename (const char *fname)
{
    char *p = strrchr(fname, '.');

    return p + 1;
}

static char *filename_from_funname (char *fname, 
				    const char *funname)
{
    build_path(fname, paths.dotdir, "pkgedit", NULL);
    strcat(fname, ".");
    strcat(fname, funname);
    return fname;
}

static void close_code_editor (const char *funname)
{
    char fname[FILENAME_MAX];
    GtkWidget *w;

    filename_from_funname(fname, funname);
    w = match_window_by_filename(fname);

    if (w != NULL) {
	gtk_widget_destroy(w);
    }    
}

static void finfo_free (function_info *finfo)
{
    int i;

    g_free(finfo->fname);
    g_free(finfo->author);
    g_free(finfo->version);
    g_free(finfo->date);
    g_free(finfo->pkgdesc);
    g_free(finfo->sample);
    g_free(finfo->help);

    if (finfo->pubnames != NULL) {
	for (i=0; i<finfo->n_pub; i++) {
	    close_code_editor(finfo->pubnames[i]);
	}
	free_strings_array(finfo->pubnames, finfo->n_pub);
    }
	
    if (finfo->privnames != NULL) {
	for (i=0; i<finfo->n_priv; i++) {
	    close_code_editor(finfo->privnames[i]);
	}
	free_strings_array(finfo->privnames, finfo->n_priv);
    }

    if (finfo->samplewin != NULL) {
	gtk_widget_destroy(finfo->samplewin->main);
    }

    if (finfo->popup != NULL) {
	gtk_widget_destroy(finfo->popup);
    }

    free(finfo);
}

static void finfo_set_modified (function_info *finfo, gboolean s)
{
    finfo->modified = s;
    gtk_widget_set_sensitive(finfo->save, s);
}

static void login_init_or_free (login_info *linfo, int freeit)
{
    static gchar *login;
    static gchar *pass;

    if (freeit) {
	if (!linfo->canceled) {
	    g_free(login);
	    g_free(pass);
	    login = g_strdup(linfo->login);
	    pass = g_strdup(linfo->pass);
	}
	g_free(linfo->login);
	g_free(linfo->pass);
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

static void login_finalize (GtkWidget *w, login_info *linfo)
{
    linfo->login = entry_box_get_trimmed_text(linfo->login_entry);
    linfo->pass = entry_box_get_trimmed_text(linfo->pass_entry);

    gtk_widget_destroy(linfo->dlg);
}

/* Used by the File, Save dialog when saving a package de novo:
   we construct a name bsed on the (first) public interface
   in the package; the user is free to change this.
*/

void get_default_package_name (char *fname, gpointer p, int mode)
{
    function_info *finfo = (function_info *) p;

    strcpy(fname, finfo->pubnames[0]);

    if (mode == SAVE_FUNCTIONS_AS) {
	strcat(fname, ".inp");
    } else {
	strcat(fname, ".gfn");
    }	
}

/* Check the user-supplied version string for the package: shoould be
   something like "1" or "1.2" or maybe "1.2.3" 
*/

static int check_version_string (const char *s)
{
    int dotcount = 0;
    int err = 0;

    /* must start and end with a digit */
    if (!isdigit(*s) || (*s && !isdigit(s[strlen(s) - 1]))) {
	err = 1;
    }

    while (*s && !err) {
	if (*s == '.' && ++dotcount > 2) {
	    /* max of two dots exceeded */
	    err = 1;
	} else if (!isdigit(*s) && strspn(s, ".") != 1) {
	    /* no more than one dot in a row allowed */
	    err = 1;
	}
	s++;
    }

    return err;
}

/* Callback from the Save button when editing a function package.  We
   assemble the relevant info then if the package is new and has not
   been saved yet (which is flagged by finfo->fname being NULL) we
   offer a file selector, otherwise we go ahead and save using the
   package's recorded filename.  
*/

static void finfo_save (GtkWidget *w, function_info *finfo)
{
    char **fields[] = {
	&finfo->author,
	&finfo->version,
	&finfo->date,
	&finfo->pkgdesc
    };
    int i, err = 0;

    for (i=0; i<NENTRIES && !err; i++) {
	g_free(*fields[i]);
	*fields[i] = entry_box_get_trimmed_text(finfo->entries[i]);
	if (*fields[i] == NULL || !strcmp(*fields[i], "missing")) {
	    gtk_entry_set_text(GTK_ENTRY(finfo->entries[i]), "missing");
	    gtk_editable_select_region(GTK_EDITABLE(finfo->entries[i]), 0, -1);
	    gtk_widget_grab_focus(finfo->entries[i]);
	    return;
	} 
    }

    g_free(finfo->help);
    finfo->help = textview_get_trimmed_text(finfo->text);
    if (finfo->help == NULL || !strcmp(finfo->help, "missing")) {
	textview_set_text_selected(finfo->text, "missing");
	gtk_widget_grab_focus(finfo->text);
	return;
    }

    if (finfo->sample == NULL) {
	infobox(_("Please add a sample script for this package"));
	return;
    }

    if (check_version_string(finfo->version)) {
	errbox(_("Invalid version string: use numbers and '.' only"));
	return;
    }

    if (finfo->fname == NULL) {
	file_selector_with_parent(SAVE_FUNCTIONS, FSEL_DATA_MISC, 
				  finfo, finfo->dlg);
    } else {
	err = save_function_package(finfo->fname, finfo);
    }
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

static gboolean update_active_func (GtkComboBox *menu, 
				    function_info *finfo)
{
    int i = gtk_combo_box_get_active(menu);

    if (i < finfo->n_pub) {
	finfo->active = finfo->pubnames[i];
    } else {
	finfo->active = finfo->privnames[i-finfo->n_pub];
    }

    return FALSE;
}

/* callback used when editing a function in the context of
   the package editor */

int update_func_code (windata_t *vwin)
{
    char *funname;
    int err = 0;

    funname = funname_from_filename(vwin->fname);
    err = update_function_from_script(funname, vwin->fname);

    if (err) {
	gui_errmsg(err);
    } else {
	function_info *finfo = vwin->data;

	finfo_set_modified(finfo, TRUE);
    }

    return err;
}

static void edit_code_callback (GtkWidget *w, function_info *finfo)
{
    char fname[FILENAME_MAX];
    ufunc *fun;
    windata_t *vwin;
    GtkWidget *orig;
    PRN *prn = NULL;

    if (finfo->active == NULL) {
	return;
    }

    filename_from_funname(fname, finfo->active);

    orig = match_window_by_filename(fname);

    if (orig != NULL) {
	gtk_window_present(GTK_WINDOW(orig));
	return;
    }

    fun = get_user_function_by_name(finfo->active);
    if (fun == NULL) {
	errbox("Can't find the function '%s'", finfo->active);
    }

    if (bufopen(&prn)) {
	return;
    }    

    gretl_function_print_code(fun, prn);

    vwin = view_buffer(prn, 78, 350, finfo->active, EDIT_FUNC_CODE, finfo);

    if (vwin != NULL) {
	strcpy(vwin->fname, fname);
    }
}

void update_sample_script (windata_t *vwin)
{
    function_info *finfo;

    finfo = g_object_get_data(G_OBJECT(vwin->main), "finfo");

    if (finfo != NULL) {
	gchar *text = textview_get_text(vwin->text);

	free(finfo->sample);
	finfo->sample = gretl_strdup(text);
	g_free(text);
	finfo_set_modified(finfo, TRUE);
    }
}

static void
nullify_sample_window (GtkWidget *w, function_info *finfo)
{
    finfo->samplewin = NULL;
}

static void edit_sample_callback (GtkWidget *w, function_info *finfo)
{
    const char *pkgname = NULL;
    gchar *title;
    PRN *prn = NULL;

    if (finfo->samplewin != NULL) {
	gtk_window_present(GTK_WINDOW(finfo->samplewin->main));
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    if (finfo->pkg != NULL) {
	/* FIXME */
	gretl_function_get_info(finfo->pubnames[0], "pkgname", &pkgname);
    } 

    if (pkgname != NULL) {
	title = g_strdup_printf("%s-sample", pkgname);
    } else {
	title = g_strdup("sample script");
    }

    if (finfo->sample != NULL) {
	pputs(prn, finfo->sample);
	pputc(prn, '\n');
    } 

    finfo->samplewin = view_buffer(prn, 78, 350, title,
				   EDIT_FUNC_CODE, finfo);

    g_object_set_data(G_OBJECT(finfo->samplewin->main), "finfo",
		      finfo);
    g_signal_connect(G_OBJECT(finfo->samplewin->main), "destroy",
		     G_CALLBACK(nullify_sample_window), finfo);

    g_free(title);
}

static void gfn_to_script_callback (GtkWidget *w, function_info *finfo)
{
    if (finfo->n_pub + finfo->n_priv == 0) {
	warnbox("No code to save");
	return;
    }

    file_selector_with_parent(SAVE_FUNCTIONS_AS, FSEL_DATA_MISC, 
			      finfo, finfo->dlg);
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

static char **get_function_names (const int *list, int *err)
{
    char **names = NULL;
    const char *funname;
    int i, n = 0;

    for (i=1; i<=list[0] && !*err; i++) {
	funname = user_function_name_by_index(list[i]);
	if (funname == NULL) {
	    *err = E_DATA;
	} else {
	    *err = strings_array_add(&names, &n, funname);
	}
    }

    if (*err) {
	free_strings_array(names, n);
	names = NULL;
    }

    return names;
}

static int finfo_add_function_names (function_info *finfo,
				     const int *publist,
				     const int *privlist)
{
    int err = 0;

    finfo->pubnames = get_function_names(publist, &err);

    if (!err && privlist != NULL) {
	finfo->privnames = get_function_names(privlist, &err);
    }

    if (!err) {
	finfo->n_pub = publist[0];
	finfo->active = finfo->pubnames[0];
	if (privlist != NULL) {
	    finfo->n_priv = privlist[0];
	}
    }

    return err;
}

static GtkWidget *active_func_selector (function_info *finfo)
{
    GtkWidget *ifmenu;
    int i;

    ifmenu = gtk_combo_box_new_text();

    for (i=0; i<finfo->n_pub; i++) {
	gtk_combo_box_append_text(GTK_COMBO_BOX(ifmenu), finfo->pubnames[i]);
    }

    for (i=0; i<finfo->n_priv; i++) {
	gchar *s = g_strdup_printf("%s (%s)", finfo->privnames[i], "private");

	gtk_combo_box_append_text(GTK_COMBO_BOX(ifmenu), s);
	g_free(s);
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(ifmenu), 0);
    gtk_widget_show_all(ifmenu);

    return ifmenu;
}

static void dreq_select (GtkComboBox *menu, function_info *finfo)
{
    finfo->dreq = gtk_combo_box_get_active(menu);
    finfo_set_modified(finfo, TRUE);
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
    GtkWidget *datamenu, *tmp;
    int j;

    tmp = gtk_label_new(_("Data requirement"));
    gtk_misc_set_alignment(GTK_MISC(tmp), 1.0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
    gtk_widget_show(tmp);

    datamenu = gtk_combo_box_new_text();
    for (j=0; j<=FN_NEEDS_PANEL; j++) {
	gtk_combo_box_append_text(GTK_COMBO_BOX(datamenu), _(datareq[j]));
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(datamenu), finfo->dreq);

    tmp = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(tmp), datamenu, FALSE, FALSE, 0);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 1, 2, i, i+1);
    gtk_widget_show_all(tmp);

    g_signal_connect(G_OBJECT(datamenu), "changed",
		     G_CALLBACK(dreq_select), finfo);
}

static void get_maj_min_pl (int v, int *maj, int *min, int *pl)
{
    *maj = v / 10000;
    *min = (v - *maj * 10000) / 100;
    *pl = v % 10;
}

static void adjust_minver (GtkWidget *w, function_info *finfo)
{
    int val = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int lev = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "level"));
    int maj, min, pl;

    get_maj_min_pl(finfo->minver, &maj, &min, &pl);

    if (lev == 1) {
	finfo->minver = 10000 * val + 100 * min + pl;
    } else if (lev == 2) {
	finfo->minver = 10000 * maj + 100 * val + pl;
    } else if (lev == 3) {
	finfo->minver = 10000 * maj + 100 * min + val;
    }

    finfo_set_modified(finfo, TRUE);
}

static void add_minver_selector (GtkWidget *tbl, int i, 
				 function_info *finfo)
{
    GtkWidget *tmp, *spin, *hbox;
    int maj, min, pl;

    get_maj_min_pl(finfo->minver, &maj, &min, &pl);

    tmp = gtk_label_new(_("Minimum gretl version"));
    gtk_misc_set_alignment(GTK_MISC(tmp), 1.0, 0.5);
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

static void pkg_changed (gpointer p, function_info *finfo)
{
    finfo_set_modified(finfo, TRUE);
}

static GtkWidget *editable_text_box (GtkTextBuffer **pbuf)
{
    GtkTextBuffer *tbuf = gretl_text_buf_new();
    GtkWidget *w = gtk_text_view_new_with_buffer(tbuf);

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(w), 4);
    gtk_widget_modify_font(GTK_WIDGET(w), fixed_font);
    gtk_text_view_set_editable(GTK_TEXT_VIEW(w), TRUE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w), TRUE);

    *pbuf = tbuf;

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

static void insert_today (GtkWidget *w, GtkWidget *entry)
{
    gtk_entry_set_text(GTK_ENTRY(entry), print_today());    
}

static gint today_popup (GtkWidget *entry, GdkEventButton *event,
			 GtkWidget **popup)
{
    GdkModifierType mods = parent_get_pointer_mask(entry);

    if (mods & GDK_BUTTON3_MASK) {
	if (*popup == NULL) {
	    GtkWidget *menu = gtk_menu_new();
	    GtkWidget *item;

	    item = gtk_menu_item_new_with_label(_("Insert today's date"));
	    g_signal_connect(G_OBJECT(item), "activate",
			     G_CALLBACK(insert_today), entry);
	    gtk_widget_show(item);
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	    *popup = menu;
	}

	gtk_menu_popup(GTK_MENU(*popup), NULL, NULL, NULL, NULL, 
		       event->button, event->time);
    }

    return TRUE;
}

static gint query_save_package (GtkWidget *w, GdkEvent *event, 
				function_info *finfo)
{
    if (finfo->modified) {
	int resp = yes_no_dialog("gretl", _("Save changes?"), 1);

	if (resp == GRETL_CANCEL) {
	    return TRUE;
	} else if (resp == GRETL_YES) {
	    finfo_save(NULL, finfo);
	}
    }

    return FALSE;
}

static void delete_pkg_editor (GtkWidget *widget, function_info *finfo) 
{
    gint resp = 0;

    if (finfo->modified) {
	resp = query_save_package(NULL, NULL, finfo);
    }

    if (!resp) {
	gtk_widget_destroy(finfo->dlg); 
    }
}

static void toggle_upload (GtkToggleButton *b, function_info *finfo)
{
    finfo->upload = gtk_toggle_button_get_active(b);
}

/* Dialog for editing function package.  The user can get here in
   either of two ways: after selecting functions to put into a
   newly created package, or upon selecting an existing package
   for editing.  In either case it's important that the "publist"
   member of the finfo struct is not NULL.
*/

static void finfo_dialog (function_info *finfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *vbox, *hbox;
    GtkTextBuffer *hbuf = NULL;
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
    int focused = 0;
    const char *hlp;
    gchar *ltxt;
    int i;

    finfo->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_default_size(GTK_WINDOW(finfo->dlg), 640, 500);

    if (finfo->fname != NULL) {
	gchar *title = title_from_filename(finfo->fname);

	gtk_window_set_title(GTK_WINDOW(finfo->dlg), title);
	g_free(title);
    } else {
	gtk_window_set_title(GTK_WINDOW(finfo->dlg), 
			     _("gretl: function package editor"));
    } 

    g_signal_connect(G_OBJECT(finfo->dlg), "delete-event",
		     G_CALLBACK(query_save_package), finfo);
    g_signal_connect(G_OBJECT(finfo->dlg), "destroy", 
		     G_CALLBACK(finfo_destroy), finfo);

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);
    gtk_container_add(GTK_CONTAINER(finfo->dlg), vbox);
    gtk_widget_show(vbox);
			 
    tbl = gtk_table_new(NENTRIES + 1, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 4);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 5);

    for (i=0; i<NENTRIES; i++) {
	GtkWidget *entry;

	label = gtk_label_new(_(entry_labels[i]));
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, i, i+1);
	gtk_widget_show(label);

	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 40);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	gtk_widget_show(entry); 

	finfo->entries[i] = entry;

	if (entry_texts[i] != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(entry), entry_texts[i]);
	    if (i == 2) {
		g_signal_connect(G_OBJECT(entry), "button-press-event",
				 G_CALLBACK(today_popup), &finfo->popup);
	    }
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

	g_signal_connect(GTK_EDITABLE(entry), "changed",
			 G_CALLBACK(pkg_changed), finfo);
    }

    add_minver_selector(tbl, i++, finfo);
    add_data_requirement_menu(tbl, i, finfo);
    gtk_widget_show(tbl);

    if (finfo->pkg != NULL) {
	const char *pkgname = function_package_get_name(finfo->pkg);

	ltxt = g_strdup_printf(_("Help text for %s:"), pkgname);
    } else {
	ltxt = g_strdup(_("Help text:"));
    }
    hbox = label_hbox(vbox, ltxt);
    gtk_widget_show(hbox);
    g_free(ltxt);

    finfo->text = editable_text_box(&hbuf);
    text_table_setup(vbox, finfo->text);

    gretl_function_get_info(finfo->pubnames[0], "help", &hlp);
    textview_set_text(finfo->text, hlp);
    g_signal_connect(G_OBJECT(hbuf), "changed", 
		     G_CALLBACK(pkg_changed), finfo);

    /* edit code button, possibly with selector */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    button = gtk_button_new_with_label(_("Edit function code"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_code_callback), finfo);

    if (finfo->n_pub + finfo->n_priv > 1) {
	finfo->codesel = active_func_selector(finfo);
	gtk_box_pack_start(GTK_BOX(hbox), finfo->codesel, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(finfo->codesel), "changed",
			 G_CALLBACK(update_active_func), finfo);
    } 

    button = gtk_button_new_with_label(_("Save as script"));
    gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(gfn_to_script_callback), finfo);

    gtk_widget_show_all(hbox);

    /* edit sample script button */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    button = gtk_button_new_with_label(_("Edit sample script"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_sample_callback), finfo);
    gtk_widget_show_all(hbox);

    /* check box for upload option */
    button = button_in_hbox(vbox, CHECK_BUTTON, 
			    _("Upload package to server on save"));
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(toggle_upload), finfo);

    /* control button area */
    hbox = gtk_hbutton_box_new();
    gtk_button_box_set_layout(GTK_BUTTON_BOX(hbox), GTK_BUTTONBOX_END);
    gtk_box_set_spacing(GTK_BOX(hbox), 10);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* Save button */
    button = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(finfo_save), finfo);
    gtk_widget_set_sensitive(button, finfo->fname == NULL);
    finfo->save = button;

    /* Close button */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(delete_pkg_editor), finfo);

    gtk_widget_show_all(hbox);

    gtk_widget_show(finfo->dlg);
}

static void web_get_login (GtkWidget *w, gpointer p)
{
    browser_open("http://gretl.ecn.wfu.edu/apply/");
}

static void login_dialog (login_info *linfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *vbox, *hbox;
    int i;

    login_init(linfo);

    linfo->dlg = gretl_dialog_new(_("gretl: upload"), NULL, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(linfo->dlg));

    hbox = label_hbox(vbox, _("Upload function package"));
    gtk_widget_show(hbox);

    tbl = gtk_table_new(2, 2, FALSE);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 5);

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

    hbox = label_hbox(vbox, 
		      _("If you don't have a login to the gretl server\n"
			"please see http://gretl.ecn.wfu.edu/apply/.\n"
			"The 'Website' button below should open this page\n"
			"in your web browser."));
    gtk_widget_show(hbox);

    /* control button area */

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(linfo->dlg));

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
    int error_printed = 0;
    int err = 0;

    login_dialog(&linfo);

    if (linfo.canceled) {
	linfo_free(&linfo);
	return;
    }

    /* set waiting cursor */
    disp = gdk_display_get_default();
    w1 = gdk_display_get_window_at_pointer(disp, &x, &y);
    if (w1 != NULL) {
	cursor = gdk_cursor_new(GDK_WATCH);
	gdk_window_set_cursor(w1, cursor);
	gdk_display_sync(disp);
	gdk_cursor_unref(cursor);
    }

    err = gretl_file_get_contents(fname, &buf);

    if (err) {
	error_printed = 1;
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
    if (w1 != NULL) {
	gdk_window_set_cursor(w1, NULL);
    }

    if (err) {
	if (!error_printed) {
	    gui_errmsg(err);
	}
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

/* the basename of a function package file must meet some sanity
   requirements */

static int check_package_filename (function_info *finfo, const char *fname)
{
    const char *p = strrchr(fname, SLASH);
    int n, err = 0;

    /* get basename in p */
    if (p == NULL) {
	p = fname;
    } else {
	p++;
    }

    if (!has_suffix(p, ".gfn")) {
	/* must have the right suffix */
	err = 1;
    } else {
	n = strlen(p) - 4;
	if (n >= FN_NAMELEN) {
	    /* too long */
	    err = 1;
	} else if (gretl_namechar_spn(p) != n) {
	    /* contains funny stuff */
	    err = 1;
	}
    }

    if (err) {
	errbox("Invalid package filename: the name must be less than\n"
	       "32 characters in length, must include only ASCII letters,\n"
	       "numbers and '_', and must end with \".gfn\"");
    }

    return err;
}

/* callback from file selector when saving a function package */

int save_function_package (const char *fname, gpointer p)
{
    function_info *finfo = p;
    int err = 0;

    if (finfo->fname == NULL) {
	/* no filename recorded yet */
	err = check_package_filename(finfo, fname);
	if (err) {
	    return err;
	}
	finfo->fname = g_strdup(fname);
    } 

    gretl_function_set_info(finfo->pubnames[0], finfo->help);

    if (finfo->pkg == NULL) {
	/* starting from scratch */
	finfo->pkg = function_package_new(fname, finfo->pubnames, finfo->n_pub,
					  finfo->privnames, finfo->n_priv,
					  &err);
    } else {
	err = function_package_connect_funcs(finfo->pkg, finfo->pubnames, finfo->n_pub,
					     finfo->privnames, finfo->n_priv);
    }

    if (!err) {
	err = function_package_set_properties(finfo->pkg,
					      "author",  finfo->author,
					      "version", finfo->version,
					      "date",    finfo->date,
					      "description", finfo->pkgdesc,
					      "sample-script", finfo->sample,
					      "data-requirement", finfo->dreq,
					      "min-version", finfo->minver,
					      NULL);
    }

    /* Note: if we allow "save as..." for existing, named function
       packages then we'll need a function to reset the package
       ID, in gretl_func.c.
    */

    if (!err) {
	err = write_function_package(finfo->pkg);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	const char *pname = function_package_get_name(finfo->pkg);
	gchar *title;

	title = g_strdup_printf("gretl: %s", pname);
	gtk_window_set_title(GTK_WINDOW(finfo->dlg), title);
	g_free(title);
	finfo_set_modified(finfo, FALSE);
	maybe_update_func_files_window(1);
	if (finfo->upload) {
	    do_upload(fname);
	}
    }

    return err;
}

/* callback from file selector when exporting a package in the form
   of a regular script */

int save_function_package_as_script (const char *fname, gpointer p)
{
    function_info *finfo = p;
    ufunc *fun;
    PRN *prn;
    int i, err = 0;

    prn = gretl_print_new_with_filename(fname, &err);
    if (err) {
	file_write_errbox(fname);
	return err;
    }

    pprintf(prn, "# author='%s'\n", finfo->author);
    pprintf(prn, "# version='%s'\n", finfo->version);
    pprintf(prn, "# date='%s'\n", finfo->date);

    for (i=0; i<finfo->n_priv; i++) {
	fun = get_user_function_by_name(finfo->privnames[i]);
	if (fun != NULL) {
	    pputc(prn, '\n');
	    gretl_function_print_code(fun, prn);
	}
    }

    for (i=0; i<finfo->n_pub; i++) {
	fun = get_user_function_by_name(finfo->pubnames[i]);
	if (fun != NULL) {
	    pputc(prn, '\n');
	    gretl_function_print_code(fun, prn);
	}
    }

    if (finfo->sample != NULL) {
	int n = strlen(finfo->sample);

	pputs(prn, "\n# sample function call\n");
	pputs(prn, finfo->sample);
	if (finfo->sample[n-1] != '\n') {
	    pputc(prn, '\n');
	}
    }

    gretl_print_destroy(prn);

    return 0;
}

/* called from function selection dialog: a set of functions has been
   selected and now we need to add info on author, version, etc.
*/

void edit_new_function_package (void)
{
    function_info *finfo;
    int *publist = NULL;
    int *privlist = NULL;
    int *list = NULL;
    int err = 0;

    if (storelist == NULL) {
	return;
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    list = gretl_list_from_string(storelist, &err);

    if (!err && (list == NULL || list[0] == 0)) {
	gretl_errmsg_set("No functions are selected");
	err = E_DATA;
    }

    if (!err) {
	if (gretl_list_has_separator(list)) {
	    err = gretl_list_split_on_separator(list, &publist,
						&privlist);
	} else {
	    publist = list;
	    list = NULL;
	}
    }

    free(list);

    if (!err) {
	err = finfo_add_function_names(finfo, publist, privlist);
	if (err) {
	    gretl_errmsg_set("Couldn't find the requested functions");
	}
    }

    free(publist);
    free(privlist);

    if (err) {
	gui_errmsg(err);
	free(finfo);
    } else {
	/* set up dialog to do the actual editing */
	finfo_dialog(finfo);
    }
}

void edit_function_package (const char *fname, int *loaderr)
{
    function_info *finfo;
    int *publist = NULL;
    int *privlist = NULL;
    const char *p;
    int err = 0;

    if (!function_package_is_loaded(fname)) {
	err = load_user_function_file(fname);
	if (err) {
	    fprintf(stderr, "load_user_function_file: failed on %s\n", fname);
	    file_read_errbox(fname);
	    if (loaderr != NULL) {
		*loaderr = 1;
	    }
	    return;
	}
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    err = function_package_get_properties(fname,
					  "package",  &finfo->pkg,
					  "publist",  &publist,
					  "privlist", &privlist,
					  "author",   &finfo->author,
					  "version",  &finfo->version,
					  "date",     &finfo->date,
					  "description", &finfo->pkgdesc,
					  "sample-script", &finfo->sample,
					  "data-requirement", &finfo->dreq,
					  "min-version", &finfo->minver,
					  NULL);

    if (!err && publist != NULL) {
	err = finfo_add_function_names(finfo, publist, privlist);
    }

    free(publist);
    free(privlist);

    if (err || publist == NULL) {
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

    finfo->fname = g_strdup(fname);

    finfo_dialog(finfo);
}

void edit_package_at_startup (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	file_read_errbox(fname);
    } else {
	fclose(fp);
	edit_function_package(fname, NULL);
    }
}

int no_user_functions_check (void)
{
    int err = 0;

    if (n_free_functions() == 0) {
	int resp;

	err = 1;
	resp = yes_no_dialog(_("gretl: function packages"),
			     _("No functions are available for packaging at present.\n"
			       "Do you want to write a function now?"),
			     0);
	if (resp == GRETL_YES) {
	    do_new_script(FUNC);
	}
    } 

    return err;
}
