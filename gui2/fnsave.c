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
#include "gretl_xml.h"
#include "libset.h"
#include "dlgutils.h"
#include "datafiles.h"
#include "textbuf.h"
#include "fileselect.h"
#include "filelists.h"
#include "gretl_www.h"
#include "gretl_zip.h"
#include "winstack.h"
#include "menustate.h"
#include "selector.h"
#include "textutil.h"
#include "fncall.h"
#include "fnsave.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#include "gretl_func.h"
#include "gretl_typemap.h"

#define PKG_DEBUG 0
#define NENTRIES 5
#define N_FILE_ENTRIES 4
#define N_SPECIALS (UFUN_ROLE_MAX - 1)

enum {
    NO_WINDOW,
    MAIN_WINDOW,
    MODEL_WINDOW
};

enum {
    APPEND_SAMPLE  = 1 << 0, /* write functions as .inp: append sample? */
    WRITE_SAMPFILE = 1 << 1, /* write .spec file: also write sample script? */
    WRITE_HELPFILE = 1 << 2, /* write .spec file: also write help file? */
    WRITE_GUI_HELP = 1 << 3  /* write .spec file: also write gui help? */
} PkgSaveFlags;

typedef struct function_info_ function_info;
typedef struct login_info_ login_info;

/* package-editing dialog and associated info */

struct function_info_ {
    GtkWidget *dlg;        /* editing dialog box */
    /* entry boxes for author, etc. */
    GtkWidget *entries[NENTRIES];
    /* entry boxes for data file names */
    GtkWidget *file_entries[N_FILE_ENTRIES];
    GtkWidget *text;       /* text box for help */
    GtkWidget *codesel;    /* code-editing selector */
    GtkWidget *popup;      /* popup menu */
    GtkWidget *extra;      /* extra properties child dialog */
    GtkWidget *maintree;   /* main menu selection tree */
    GtkWidget *modeltree;  /* model menu selection tree */
    GtkWidget *currtree;   /* currently displayed menu treeview */
    GtkWidget *alttree;    /* currently undisplayed menu treeview */
    GtkWidget *treewin;    /* scrolled window to hold menu trees */
    GtkWidget *specdlg;    /* pkg spec save dialog */
    GtkWidget *validate;   /* "Validate" gfn button */
    windata_t *samplewin;  /* window for editing sample script */
    windata_t *gui_helpwin;  /* window for editing GUI help text */
    GtkUIManager *ui;      /* for dialog File menu */
    GList *codewins;       /* list of windows editing function code */
    fnpkg *pkg;            /* pointer to package being edited */
    gchar *pkgname;        /* suggested savename (temporary) */
    gchar *fname;          /* package filename */
    gchar *author;         /* package author */
    gchar *email;          /* author's email address */
    gchar *version;        /* package version number */
    gchar *date;           /* package last-revised date */
    gchar *pkgdesc;        /* package description */
    gchar *sample;         /* sample script for package */
    gchar *help;           /* package help text */
    gchar *gui_help;       /* GUI-specific help text */
    gchar *sample_fname;   /* filename: sample script */
    gchar *help_fname;     /* filename: help text */
    gchar *gui_help_fname; /* filename: GUI-specific help */
    char **pubnames;       /* names of public functions */
    char **privnames;      /* names of private functions */
    char **specials;       /* names of special functions */
    char **datafiles;      /* names of included data files */
    int n_pub;             /* number of public functions */
    int n_priv;            /* number of private functions */
    int n_files;           /* number of include data files */
    gboolean uses_subdir;  /* the package has its own subdir (0/1) */
    gchar *menupath;       /* path for menu attachment, if any */
    gchar *menulabel;      /* label for menu attachment, if any */
    int menuwin;           /* code for none/main/model window */
    char *active;          /* name of 'active' function */
    FuncDataReq dreq;      /* data requirement of package */
    int minver;            /* minimum gretl version, package */
    gboolean modified;     /* anything changed in package? */
    int save_flags;        /* see PkgSaveFlags */
    unsigned char gui_attrs[N_SPECIALS]; /* attribute flags for special funcs */
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

static int validate_package_file (const char *fname,
				  int verbose);
static void finfo_set_menuwin (function_info *finfo);
static gint query_save_package (GtkWidget *w, GdkEvent *event, 
				function_info *finfo);
static int finfo_save (function_info *finfo);
static void gfn_to_script_callback (function_info *finfo);
static void gfn_to_spec_callback (function_info *finfo);
static void do_pkg_upload (function_info *finfo);
static int pkg_has_pdf_doc (function_info *finfo);
static void edit_code_callback (GtkWidget *w, function_info *finfo);
static int check_package_filename (const char *fname,
				   int fullpath,
				   GtkWidget *parent);

function_info *finfo_new (void)
{
    function_info *finfo;

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) {
	return NULL;
    }

    finfo->specials = strings_array_new(N_SPECIALS);
    if (finfo->specials == NULL) {
	free(finfo);
	return NULL;
    }

    memset(finfo->gui_attrs, 0, N_SPECIALS);

    finfo->pkg = NULL;
    finfo->fname = NULL;
    finfo->pkgname = NULL;
    finfo->author = NULL;
    finfo->email = NULL;
    finfo->version = NULL;
    finfo->date = NULL;
    finfo->pkgdesc = NULL;
    finfo->sample = NULL;
    finfo->menupath = NULL;
    finfo->menulabel = NULL;
    finfo->menuwin = 0;

    finfo->currtree = NULL;
    finfo->alttree = NULL;
    finfo->treewin = NULL;
    finfo->ui = NULL;

    finfo->modified = FALSE;
    finfo->save_flags = WRITE_SAMPFILE |
	WRITE_HELPFILE | WRITE_GUI_HELP;

    finfo->active = NULL;
    finfo->samplewin = NULL;
    finfo->gui_helpwin = NULL;
    finfo->codewins = NULL;
    finfo->codesel = NULL;
    finfo->popup = NULL;
    finfo->extra = NULL;
    finfo->specdlg = NULL;

    finfo->help = NULL;
    finfo->gui_help = NULL;

    finfo->sample_fname = NULL;
    finfo->help_fname = NULL;
    finfo->gui_help_fname = NULL;
    
    finfo->pubnames = NULL;
    finfo->privnames = NULL;
    finfo->datafiles = NULL;

    finfo->n_pub = 0;
    finfo->n_priv = 0;
    finfo->n_files = 0;

    finfo->dreq = 0;
    finfo->minver = 10900;
    finfo->uses_subdir = 0;

    return finfo;
}

static const char *funname_from_filename (const char *fname)
{
    const char *p = strrchr(fname, '.');

    return p + 1;
}

static char *filename_from_funname (char *fname, 
				    const char *funname)
{
    build_path(fname, gretl_dotdir(), "pkgedit", NULL);
    strcat(fname, ".");
    strcat(fname, funname);
    return fname;
}

static void destroy_code_window (windata_t *vwin)
{
    gtk_widget_destroy(vwin->main);
}

static void finfo_free (function_info *finfo)
{
    g_free(finfo->fname);
    g_free(finfo->author);
    g_free(finfo->email);
    g_free(finfo->version);
    g_free(finfo->date);
    g_free(finfo->pkgdesc);
    g_free(finfo->sample);
    g_free(finfo->help);
    g_free(finfo->gui_help);

    g_free(finfo->sample_fname);
    g_free(finfo->help_fname);
    g_free(finfo->gui_help_fname);
    
    g_free(finfo->menupath);
    g_free(finfo->menulabel);

    if (finfo->pubnames != NULL) {
	strings_array_free(finfo->pubnames, finfo->n_pub);
    }
	
    if (finfo->privnames != NULL) {
	strings_array_free(finfo->privnames, finfo->n_priv);
    }

    if (finfo->specials != NULL) {
	strings_array_free(finfo->specials, N_SPECIALS);
    }

    if (finfo->datafiles != NULL) {
	strings_array_free(finfo->datafiles, finfo->n_files);
    }    

    if (finfo->samplewin != NULL) {
	gtk_widget_destroy(finfo->samplewin->main);
    }

    if (finfo->gui_helpwin != NULL) {
	gtk_widget_destroy(finfo->gui_helpwin->main);
    }    

    if (finfo->codewins != NULL) {
	g_list_foreach(finfo->codewins, (GFunc) destroy_code_window, NULL);
	g_list_free(finfo->codewins);
    }

    if (finfo->ui != NULL) {
	g_object_unref(finfo->ui);
    }

    if (finfo->popup != NULL) {
	gtk_widget_destroy(finfo->popup);
    }

    free(finfo);
}

static void pkg_save_action (GtkAction *action, function_info *finfo)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "Save")) {
	finfo_save(finfo);
    } else if (!strcmp(s, "SaveZip")) {
	file_selector_with_parent(SAVE_GFN_ZIP, FSEL_DATA_MISC,
				  finfo, finfo->dlg);
    } else if (!strcmp(s, "WriteInp")) {
	gfn_to_script_callback(finfo);
    } else if (!strcmp(s, "WriteSpec")) {
	gfn_to_spec_callback(finfo);
    } else if (!strcmp(s, "Upload")) {
	do_pkg_upload(finfo);
    }
}

const gchar *pkgsave_ui = 
    "<ui>"
    "  <popup>"
    "    <menuitem action='Save'/>"
    "    <menuitem action='SaveZip'/>"
    "    <menuitem action='WriteInp'/>"
    "    <menuitem action='WriteSpec'/>"
    "    <menuitem action='Upload'/>"
    "  </popup>"
    "</ui>";

static GtkActionEntry pkgsave_items[] = {
    { "Save", NULL, N_("_Save gfn"), NULL, NULL, G_CALLBACK(pkg_save_action) },
    { "SaveZip", NULL, N_("Save _zip file..."), NULL, NULL, G_CALLBACK(pkg_save_action) },
    { "WriteInp", NULL, N_("Save as _script..."), NULL, NULL, G_CALLBACK(pkg_save_action) },
    { "WriteSpec", NULL, N_("_Write spec file..."), NULL, NULL, G_CALLBACK(pkg_save_action) },
    { "Upload", NULL, N_("_Upload to server..."), NULL, NULL, G_CALLBACK(pkg_save_action) },
};

static void save_popup_pos (GtkMenu *menu,
			    gint *x,
			    gint *y,
			    gboolean *push_in,
			    gpointer data)
{
    GtkWidget *button = data;
    gint wx, wy, tx, ty;

    gdk_window_get_origin(gtk_widget_get_window(button), &wx, &wy);
    gtk_widget_translate_coordinates(button, gtk_widget_get_toplevel(button), 
				     0, 0, &tx, &ty);
    *x = wx + tx - 80;
    *y = wy + ty - 128;
    *push_in = TRUE;
}

static void pkg_save_popup (GtkWidget *button,
			    function_info *finfo)
{
    GtkWidget *menu;
    gboolean cond;

    if (finfo->ui == NULL) {
	GtkActionGroup *actions;
	
	finfo->ui = gtk_ui_manager_new();
	actions = gtk_action_group_new("PkgActions");
	gtk_action_group_set_translation_domain(actions, "gretl");
	gtk_action_group_add_actions(actions, pkgsave_items,
				     G_N_ELEMENTS(pkgsave_items),
				     finfo);
	gtk_ui_manager_add_ui_from_string(finfo->ui, pkgsave_ui, -1, NULL);
	gtk_ui_manager_insert_action_group(finfo->ui, actions, 0);
	g_object_unref(actions);
    }

    /* set menu item sensitivities */
    flip(finfo->ui, "/popup/Save", finfo->modified);
    cond = finfo->fname != NULL;
    flip(finfo->ui, "/popup/Upload", cond);
    cond = pkg_has_pdf_doc(finfo) || finfo->datafiles != NULL;
    flip(finfo->ui, "/popup/SaveZip", cond);    

    menu = gtk_ui_manager_get_widget(finfo->ui, "/popup");

    gtk_menu_popup(GTK_MENU(menu), NULL, NULL,
		   save_popup_pos,
		   button, 0,
		   gtk_get_current_event_time());    
}

static void finfo_set_modified (function_info *finfo, gboolean s)
{
    finfo->modified = s;
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
	linfo->canceled = 1;
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
    if (linfo->login == NULL) {
	gtk_widget_grab_focus(linfo->login_entry);
	return;
    }

    linfo->pass = entry_box_get_trimmed_text(linfo->pass_entry);
    if (linfo->pass == NULL) {
	gtk_widget_grab_focus(linfo->pass_entry);
	g_free(linfo->login);
	return;
    }

    gtk_widget_destroy(linfo->dlg);
}

/* Used by the File, Save dialog when saving a package,
   or saving packaged functions as a script, or writing
   a .spec file based on a package.
*/

void get_default_package_name (char *fname, gpointer p, int mode)
{
    function_info *finfo = (function_info *) p;
    const char *pkgname;

    if (finfo->pkgname != NULL) {
	/* first=time save */
	pkgname = finfo->pkgname;
    } else {
	pkgname = function_package_get_name(finfo->pkg);
    }

    if (pkgname != NULL && *pkgname != '\0') {
	strcpy(fname, pkgname);
    } else if (finfo->n_pub > 0) {
	int i;

	/* We'll prefer, as definitive of the package name, the
	   first public interface name that doesn't include an
	   underscore.
	*/
	*fname = '\0';
	for (i=0; i<finfo->n_pub; i++) {
	    if (strchr(finfo->pubnames[i], '_') == NULL) {
		strcpy(fname, finfo->pubnames[i]);
		break;
	    }
	}
	/* fallback */
	if (*fname == '\0') {
	    strcpy(fname, finfo->pubnames[0]);
	}
    } else {
	strcpy(fname, "pkg");
    }

    if (mode == SAVE_FUNCTIONS_AS) {
	strcat(fname, ".inp");
    } else if (mode == SAVE_GFN_SPEC) {
	strcat(fname, ".spec");
    } else if (mode == SAVE_GFN_ZIP) {
	strcat(fname, ".zip");
    } else {
	strcat(fname, ".gfn");
    }	
}

/* fairly minimal check here! */

static int check_email_string (const char *s)
{
    int err = 0;
    
    if (strchr(s, ' ') != NULL) {
	/* no spaces allowed */
	err = 1;
    } else if (strchr(s, '@') == NULL) {
	/* must include "at"-sign */
	err = 1;
    }

    return err;
}

/* Check the user-supplied version string for the package: should be
   something like "1" or "1.02"
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
	if (!isdigit(*s) && *s != '.') {
	    /* only dots and digits allowed */
	    err = 1;
	} else if (*s == '.' && ++dotcount > 1) {
	    /* max of one dot exceeded */
	    err = 1;
	}
	s++;
    }

    return err;
}

static int pkg_has_pdf_doc (function_info *finfo)
{
    if (finfo->help != NULL &&
	!strncmp(finfo->help, "pdfdoc", 6)) {
	return 1;
    } else {
	return 0;
    }
}

static int pkg_path_is_toplevel (function_info *finfo,
				 const char *pkgname)
{
    gchar *test;
    int ret;

    /* look for pattern "functions/mypkg.gfn" */
    test = g_strdup_printf("functions%c%s.gfn", SLASH, pkgname);
    ret = strstr(finfo->fname, test) != NULL;
    g_free(test);

    return ret;
}

static void set_gfn_save_opt (GtkWidget *w, int *opt)
{
    *opt = widget_get_int(w, "action");
}

gint compare_pubnames (gconstpointer a, gconstpointer b)
{
    const char *sa = a, *sb = b;

    if (strchr(sa, '_') == NULL && strchr(sb, '_') != NULL) {
	return -1;
    } else if (strchr(sb, '_') == NULL && strchr(sa, '_') != NULL) {
	return 1;
    } else {
	return 0;
    }
}

static GList *suitable_savenames (function_info *finfo, gchar **pkgname)
{
    GList *list = NULL;
    char *s;
    int i, role;

    for (i=0; i<finfo->n_pub; i++) {
	role = i + 1;
	s = finfo->pubnames[i];
	role = strings_array_position(finfo->specials, N_SPECIALS, s) + 1;
	if (role == 0 || role == UFUN_GUI_MAIN) {
	    list = g_list_append(list, finfo->pubnames[i]);
	}
    }

    if (list != NULL) {
	list = g_list_sort(list, compare_pubnames);
	*pkgname = g_strdup((char *) list->data);
    }

    return list;
}

static void save_gfn_ok (GtkButton *button, gchar **pkgname)
{
    GtkWidget *dialog, *combo, *entry;
    const gchar *s;
    int err = 0;

    dialog = g_object_get_data(G_OBJECT(button), "dialog");
    combo = g_object_get_data(G_OBJECT(button), "combo");
    entry = gtk_bin_get_child(GTK_BIN(combo));
    s = gtk_entry_get_text(GTK_ENTRY(entry));

    if (string_is_blank(s)) {
	err = 1;
    } else {
	g_free(*pkgname);
	*pkgname = g_strdup(s);
	g_strstrip(*pkgname);
	err = check_package_filename(*pkgname, 0, dialog);
    }

    if (err) {
	gtk_widget_grab_focus(entry);
    } else {
	gtk_widget_destroy(dialog);
    }
}

static void save_gfn_cancel (GtkButton *button, int *retval)
{
    GtkWidget *dialog;

    dialog = g_object_get_data(G_OBJECT(button), "dialog");
    *retval = GRETL_CANCEL;
    gtk_widget_destroy(dialog);
}

static void save_gfn_delete (GtkWidget *w, GdkEvent *event, int *retval)
{
    *retval = GRETL_CANCEL;
}

int save_gfn_dialog (function_info *finfo, gchar **pkgname)
{
    const char *opts[] = {
	N_("Save the file to its standard \"installed\" location"),
	N_("Save it to a location of your own choosing")
    };    
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *w;
    GtkWidget *combo, *entry;
    GtkWidget *button = NULL;
    GSList *group = NULL;
    GList *flist = NULL;
    int i, ret = 0;

    if (maybe_raise_dialog()) {
	return ret;
    }

    dialog = gretl_dialog_new(NULL, finfo->dlg, GRETL_DLG_BLOCK);
    g_signal_connect(G_OBJECT(dialog), "delete-event",
		     G_CALLBACK(save_gfn_delete), &ret);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);
    w = gtk_label_new(_("Save gfn file"));
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);

    for (i=0; i<2; i++) {
	button = gtk_radio_button_new_with_label(group, _(opts[i]));
	gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
	g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(i));
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_gfn_save_opt), &ret);
	if (i == 0) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
	}
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    /* label + combo box with entry for naming the package */

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Select or type a name for the package:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    combo = combo_box_text_new_with_entry();
    flist = suitable_savenames(finfo, pkgname);
    if (flist != NULL) {
	set_combo_box_strings_from_list(combo, flist);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    }
    gtk_box_pack_end(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" */
    w = cancel_button(hbox);
    gtk_widget_set_can_default(w, FALSE);
    g_object_set_data(G_OBJECT(w), "dialog", dialog);
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(save_gfn_cancel), &ret);

    /* "OK" */
    w = ok_button(hbox);
    g_object_set_data(G_OBJECT(w), "combo", combo);
    g_object_set_data(G_OBJECT(w), "dialog", dialog);
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(save_gfn_ok), pkgname);
    gtk_widget_grab_default(w);

    entry = gtk_bin_get_child(GTK_BIN(combo));
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    gtk_widget_show_all(dialog);

    return ret;
}

static int overwrite_gfn_check (const char *fname,
				GtkWidget *parent,
				int *notified)
{
    struct stat sbuf = {0};
    int resp = GRETL_YES;

    if (gretl_stat(fname, &sbuf) == 0) {
	gchar *msg;
	
	msg = g_strdup_printf("%s\n\n%s\n%s", fname,
			      _("A file of this name already exists."),
			      _("OK to overwite it?"));
	resp = yes_no_dialog(NULL, msg, parent);
	g_free(msg);
	if (notified != NULL) {
	    *notified = 1;
	}	
    }

    return resp;
}

/* In saving a new package, the user has chosen to write to
   the "install" location. We need to figure out the correct
   path (possibly creating a package-specific directory) then
   save the gfn file. We should give the user some feedback
   whether or not this is successful.
*/

static int install_gfn (function_info *finfo, gchar *pkgname)
{
    char savepath[FILENAME_MAX];
    gchar *msg;
    int notified = 0;
    int err = 0;

    get_default_dir_for_action(savepath, SAVE_FUNCTIONS);
    
    if (finfo->uses_subdir) {
	strcat(savepath, pkgname);
	err = gretl_mkdir(savepath);
	if (err) {
	    gui_errmsg(err);
	} else {
	    slash_terminate(savepath);
	}
    }

    if (!err) {
	int resp;
	
	strcat(savepath, pkgname);
	strcat(savepath, ".gfn");
	resp = overwrite_gfn_check(savepath, finfo->dlg,
				   &notified);
	if (resp != GRETL_YES) {
	    return 0;
	}
	err = save_function_package(savepath, finfo);
    }

    if (!err && !notified) {
	msg = g_strdup_printf(_("Wrote gfn file as\n%s"), savepath);
	msgbox(msg, GTK_MESSAGE_INFO, finfo->dlg);
	g_free(msg);
    }
    
    return err;
}

static int check_help_text (function_info *finfo)
{
    const char *missing = "???";
    char *s = finfo->help;
    int err = 0;
    
    if (s == NULL || !strcmp(s, missing)) {
	warnbox(_("Please complete all fields"));
	textview_set_text_selected(finfo->text, missing);
	err = 1;
    } else if (!strncmp(s, "pdfdoc:", 7) &&
	       strstr(s + 7, ".pdf") == NULL) {
	warnbox(_("Please complete all fields"));
	err = 1;
    } else {
	/* We can't do any kind of semantic test on the help
	   text here, but we can at least check it for overly
	   long lines, which will look ugly and broken in the
	   gretl GUI.
	*/
	int ulen, n = 0, err = 0;
	char *p = s;

	while (*s && !err) {
	    if (*s == '\n') {
		p = s + 1; /* start of line marker */
		n = 0;
	    } else {
		n++;
	    }
	    if (n > 72) {
		ulen = g_utf8_strlen(p, n);
		if (ulen > 72) {
		    warnbox(_("Please limit help text lines to 72 characters"));
		    err = 1;
		}
	    }
	    s++;
	}
    }

    return err;
}

/* Callback from "Save" , when editing a function package. We first
   assemble and check the relevant info then if the package is new and
   has not been saved yet (which is flagged by finfo->fname being
   NULL) we offer a file selector, otherwise we go ahead and save
   using the package's recorded filename.
*/

static int finfo_save (function_info *finfo)
{
    const char *missing = "???";
    char **fields[] = {
	&finfo->author,
	&finfo->email,
	&finfo->version,
	&finfo->date,
	&finfo->pkgdesc
    };
    int i, err = 0;

    for (i=0; i<NENTRIES && !err; i++) {
	g_free(*fields[i]);
	*fields[i] = entry_box_get_trimmed_text(finfo->entries[i]);
	if (*fields[i] == NULL || !strcmp(*fields[i], missing)) {
	    warnbox(_("Please complete all fields"));
	    gtk_entry_set_text(GTK_ENTRY(finfo->entries[i]), missing);
	    gtk_editable_select_region(GTK_EDITABLE(finfo->entries[i]), 0, -1);
	    gtk_widget_grab_focus(finfo->entries[i]);
	    return 1;
	} 
    }

    g_free(finfo->help);
    finfo->help = textview_get_trimmed_text(finfo->text);
    
    if (check_help_text(finfo)) {
	gtk_widget_grab_focus(finfo->text);
	return 1;
    }

    if (finfo->sample == NULL) {
	warnbox(_("Please add a sample script for this package"));
	return 1;
    }

    if (check_email_string(finfo->email)) {
	errbox(_("Please supply a valid email address"));
	return 1;
    } else {
	set_author_mail(finfo->email);
    }

    if (check_version_string(finfo->version)) {
	errbox(_("Invalid version string: use numbers and '.' only"));
	return 1;
    }

    if (!finfo->uses_subdir) {
	if (finfo->n_files > 0 || pkg_has_pdf_doc(finfo)) {
	    finfo->uses_subdir = 1;
	}
    }

    if (finfo->fname == NULL) {
	/* first-time save */
	gchar *pkgname = NULL;
	int resp;

	resp = save_gfn_dialog(finfo, &pkgname);

	if (resp == 0) {
	    /* "install" the gfn file */
	    err = install_gfn(finfo, pkgname);
	} else if (resp == 1) {
	    /* we'll save to a user-determined location: but first 
	       set the default package name on @finfo
	    */
	    finfo->pkgname = pkgname;
	    pkgname = NULL;
	}
	g_free(pkgname);
	if (resp < 1) {
	    /* cancel or "install": handled above */
	    return err;
	}
    }

    if (finfo->fname == NULL) {
	/* note: the callback from the file selector is 
	   save_function_package()
	*/
	file_selector_with_parent(SAVE_FUNCTIONS, FSEL_DATA_MISC, 
				  finfo, finfo->dlg);
    } else {
	err = save_function_package(finfo->fname, finfo);
    }

    g_free(finfo->pkgname);
    finfo->pkgname = NULL;

    return err;
}

static void finfo_destroy (GtkWidget *w, function_info *finfo)
{
    finfo_free(finfo);
}

static gboolean update_active_func (GtkComboBox *menu, 
				    function_info *finfo)
{
    int i = 0;

    if (menu != NULL) {
	i = gtk_combo_box_get_active(menu);
	if (i < 0) {
	    i = 0;
	}
    }

    if (i < finfo->n_pub) {
	finfo->active = finfo->pubnames[i];
    } else {
	finfo->active = finfo->privnames[i-finfo->n_pub];
    }

    return FALSE;
}

/* Given a line "function ..." get the function name, with
   some error checking.  The @s we are given here is at an 
   offset of 9 bytes into the line, skipping "function ".
*/

static int extract_funcname (const char *s, const char *origname)
{
    char newname[FN_NAMELEN];
    char word[FN_NAMELEN];
    int n, type, err = 0;

    s += strspn(s, " ");
    n = strcspn(s, " (");

    if (n == 0 || n > FN_NAMELEN - 1) {
	return E_DATA;
    }

    *newname = *word = '\0';
    strncat(word, s, n);

    if (!strcmp(word, "void")) {
	type = GRETL_TYPE_VOID;
    } else {
	type = gretl_type_from_string(word);
    }

    if (!ok_function_return_type(type)) {
	gretl_errmsg_sprintf("%s: bad or missing return type", origname);
	err = E_DATA;
    } else {
	s += n;
	s += strspn(s, " ");
	n = strcspn(s, " (");
	if (n == 0 || n > FN_NAMELEN - 1) {
	    err = E_DATA;
	} else {
	    strncat(newname, s, n);
	    if (strcmp(newname, origname)) {
		gretl_errmsg_set(_("You can't change the name of a function here"));
		err = E_DATA;
	    }
	}
    }

    return err;
}

static int pretest_funcname (char *buf, const char *origname)
{
    char *s, line[MAXLINE];
    int err = 0;

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && !err) {
	s = line + strspn(line, " \t");
	if (!strncmp(s, "function ", 9)) {
	    err = extract_funcname(s + 9, origname);
	    break;
	}
    }

    bufgets_finalize(buf);

    return err;
}

/* callback used when editing a function in the context of the package
   editor: save window-content to file and pass this to gretl_func to
   revise the function definition.  
*/

int update_func_code (windata_t *vwin)
{
    gchar *text = textview_get_text(vwin->text);
    function_info *finfo = vwin->data;
    const char *funname;
    int err;

    funname = funname_from_filename(vwin->fname);
    err = pretest_funcname(text, funname);

    if (!err) {
	int save_batch = gretl_in_batch_mode();

	set_current_function_package(finfo->pkg);
	err = execute_script(NULL, text, NULL, INCLUDE_EXEC, NULL);
	set_current_function_package(NULL);
	gretl_set_batch_mode(save_batch);
    }

    g_free(text);

    if (err) {
	gui_errmsg(err);
    } else {
	mark_vwin_content_saved(vwin);
	finfo_set_modified(finfo, TRUE);
    }    

    return err;
}

static void finfo_remove_codewin (GtkWidget *w, function_info *finfo)
{
    gpointer p = g_object_get_data(G_OBJECT(w), "vwin");

    finfo->codewins = g_list_remove(finfo->codewins, p);
}

static gboolean funcname_limit (gunichar c, gpointer p)
{
    return (!isalpha(c) && c != '_');
}

static gint catch_codewin_key (GtkWidget *w, GdkEventKey *event, 
			       function_info *finfo)
{
    if (finfo->n_pub + finfo->n_priv < 2) {
	/* we don't have multiple functions */
	return FALSE;
    }

    /* implement Alt-dot in a given function editing window to
       traverse to another window editing a different function
    */

    if ((event->state & GDK_MOD1_MASK) && event->keyval == GDK_period) {
	/* Alt + dot */
	windata_t *vwin = g_object_get_data(G_OBJECT(w), "vwin");
	GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
	GtkTextBuffer *buf = gtk_text_view_get_buffer(view); 
	GtkTextMark *mark = gtk_text_buffer_get_insert(buf);
	GtkTextIter iter, istart, iend;
	gchar *word = NULL;

	gtk_text_buffer_get_iter_at_mark(buf, &iter, mark);
	istart = iend = iter;
	if (gtk_text_iter_backward_find_char(&istart, funcname_limit, NULL, NULL) &&
	    gtk_text_iter_forward_char(&istart) &&
	    gtk_text_iter_forward_find_char(&iend, funcname_limit, NULL, NULL)) {
	    word = gtk_text_buffer_get_text(buf, &istart, &iend, FALSE);
	}

	if (word != NULL) {
	    /* we got a "word": is it the name of a function in
	       this pachage? */
	    char *active = NULL;
	    int i;

	    for (i=0; i<finfo->n_pub && !active; i++) {
		if (!strcmp(word, finfo->pubnames[i])) {
		    active = finfo->pubnames[i];
		}
	    }
	    for (i=0; i<finfo->n_priv && !active; i++) {
		if (!strcmp(word, finfo->privnames[i])) {
		    active = finfo->privnames[i];
		}
	    }
	    if (active != NULL) {
		finfo->active = active;
		edit_code_callback(NULL, finfo);
	    }
	    g_free(word);
	}

	return TRUE;
    }

    return FALSE;
}

static void finfo_add_codewin (function_info *finfo, windata_t *vwin)
{
    finfo->codewins = g_list_append(finfo->codewins, vwin);

    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_codewin_key), finfo);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(finfo_remove_codewin),
		     finfo);
}

static windata_t *get_codewin_by_filename (const char *fname,
					   function_info *finfo)
{
    GList *list = finfo->codewins;
    windata_t *vwin;

    while (list) {
	vwin = list->data;
	if (vwin != NULL && !strcmp(fname, vwin->fname)) {
	    return vwin;
	}
	list = g_list_next(list);
    }

    return NULL;
}

/* editing a public interface or private function belonging
   to a package: callback from "Edit function code" button.
*/

static void edit_code_callback (GtkWidget *w, function_info *finfo)
{
    char *funname = finfo->active;
    char fname[FILENAME_MAX];
    ufunc *fun;
    windata_t *vwin;
    PRN *prn = NULL;

    if (funname == NULL) {
	return;
    }

    filename_from_funname(fname, funname);

    vwin = get_codewin_by_filename(fname, finfo);
    if (vwin != NULL) {
	gtk_window_present(GTK_WINDOW(vwin->main));
	return;
    }

    fun = get_function_from_package(funname, finfo->pkg);
    if (fun == NULL) {
	/* the package may not be saved yet */
	fun = get_user_function_by_name(funname);
	if (fun == NULL) {
	    errbox_printf(_("Can't find the function '%s'"), funname);
	}
    }

    if (bufopen(&prn)) {
	return;
    }    

    gretl_function_print_code(fun, tabwidth, prn);

    vwin = view_buffer(prn, SCRIPT_WIDTH, SCRIPT_HEIGHT, 
		       finfo->active, EDIT_PKG_CODE, finfo);

    if (vwin != NULL) {
	strcpy(vwin->fname, fname);
	finfo_add_codewin(finfo, vwin);
	set_window_delete_filename(vwin);
    }
}

/* used by callback from Exec in sample script editor window */

gchar *package_sample_get_script (windata_t *vwin)
{
    function_info *finfo = vwin->data;
    gchar *buf = textview_get_text(vwin->text);
    const char *pkgname;
    gchar *ret;
    char *p, line[MAXLINE];
    gsize retsize;
    int n, done = 0;

    if (buf == NULL || *buf == '\0' || finfo->pkg == NULL) {
	return buf;
    }

    pkgname = function_package_get_name(finfo->pkg);

    /* allow for adding "# ", and possibly appending a newline */
    retsize = strlen(buf) + 5;
    ret = g_malloc(retsize);
    *ret = '\0';

    /* We need to comment out "include <self>.gfn" if such a line is
       included in the sample script: the package is already in memory
       and re-loading it now may be disruptive.
    */

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	if (!done) {
	    p = line + strspn(line, " ");
	    if (!strncmp(p, "include ", 8)) {
		p += 8;
		p += strspn(p, " ");
		n = gretl_namechar_spn(p);
		if (!strncmp(p, pkgname, n)) {
		    g_strlcat(ret, "# ", retsize);
		    done = 1;
		} 
	    }
	}
	g_strlcat(ret, line, retsize);
    }

    bufgets_finalize(buf);   
    g_free(buf);

    n = strlen(ret);
    if (ret[n-1] != '\n') {
	g_strlcat(ret, "\n", retsize);
    }    

    return ret;
}

/* callback from Save in sample script editor window */

void update_sample_script (windata_t *vwin)
{
    function_info *finfo;

    finfo = g_object_get_data(G_OBJECT(vwin->main), "finfo");

    if (finfo != NULL) {
	gchar *text = textview_get_text(vwin->text);

	free(finfo->sample);
	if (text == NULL || string_is_blank(text)) {
	    finfo->sample = NULL;
	} else {
	    finfo->sample = gretl_strdup(text);
	}
	g_free(text);
	mark_vwin_content_saved(vwin);
	finfo_set_modified(finfo, TRUE);
    }
}

/* callback from Save in GUI help text editor window */

void update_gfn_gui_help (windata_t *vwin)
{
    function_info *finfo;

    finfo = g_object_get_data(G_OBJECT(vwin->main), "finfo");

    if (finfo != NULL) {
	gchar *text = textview_get_text(vwin->text);

	free(finfo->gui_help);
	if (text == NULL || string_is_blank(text)) {
	    finfo->gui_help = NULL;
	} else {
	    finfo->gui_help = gretl_strdup(text);
	}
	g_free(text);
	mark_vwin_content_saved(vwin);
	finfo_set_modified(finfo, TRUE);
    }
}

static void nullify_sample_window (GtkWidget *w, function_info *finfo)
{
    finfo->samplewin = NULL;
}

static void nullify_gui_helpwin (GtkWidget *w, function_info *finfo)
{
    finfo->gui_helpwin = NULL;
}

/* edit the sample script for a package: callback from
   "Edit sample script" button in packager
*/

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
	pkgname = function_package_get_name(finfo->pkg);
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
				   EDIT_PKG_SAMPLE, finfo);

    g_object_set_data(G_OBJECT(finfo->samplewin->main), "finfo",
		      finfo);
    g_signal_connect(G_OBJECT(finfo->samplewin->main), "destroy",
		     G_CALLBACK(nullify_sample_window), finfo);

    g_free(title);
}

/* Callback to launch dialog for adding or removing functions.
   We need to be careful here: if the package's "extra
   properties" dialog is open, its content is liable to be
   out-dated by changes in the public and/or private
   function lists. Since it would be very complicated and
   error-prone to adjust this content on the fly, we'll
   insist that the user closes the extra props dialog
   first.
*/

static void add_remove_callback (GtkWidget *w, function_info *finfo)
{
    if (finfo->extra != NULL) {
	const char *msg = N_("Before adding or removing functions, please close\n"
			     "the \"extra properties\" dialog (after applying any\n"
			     "changes you wish to keep).");

	gtk_window_present(GTK_WINDOW(finfo->extra));
	msgbox(_(msg), GTK_MESSAGE_INFO, finfo->extra);
	return;
    }

    add_remove_functions_dialog(finfo->pubnames, finfo->n_pub,
				finfo->privnames, finfo->n_priv,
				finfo->pkg, finfo);
}

static void gfn_to_script_callback (function_info *finfo)
{
    gint resp;
    
    if (finfo->n_pub + finfo->n_priv == 0) {
	warnbox("No code to save");
	return;
    }

    if (finfo->sample != NULL) {
	resp = yes_no_cancel_dialog("gretl",
				    _("Saving packaged functions as script:\n"
				      "include the sample script?"),
				    finfo->dlg);
	
	if (canceled(resp)) {
	    return;
	}

	if (resp == GRETL_YES) {
	    finfo->save_flags |= APPEND_SAMPLE;
	} else {
	    finfo->save_flags &= ~APPEND_SAMPLE;
	}
    }

    file_selector_with_parent(SAVE_FUNCTIONS_AS, FSEL_DATA_MISC, 
			      finfo, finfo->dlg);
}

struct spec_info {
    GtkWidget *dialog;
    GtkWidget *checks[3];
    GtkWidget *entries[3];
    int *flags;
    function_info *finfo;
    int retval;
};

static void reset_finfo_filename (function_info *finfo, int i, gchar *src)
{
    if (i == 0) {
	g_free(finfo->sample_fname);
	finfo->sample_fname = src;
    } else if (i == 1) {
	g_free(finfo->help_fname);
	finfo->help_fname = src;
    } else {
	g_free(finfo->gui_help_fname);
	finfo->gui_help_fname = src;
    }
}

static void spec_save_ok (GtkWidget *button, gpointer data)
{
    struct spec_info *sinfo = data;
    function_info *finfo = sinfo->finfo;
    gchar *fname;
    int i, flag;

    for (i=0; i<3; i++) {
	if (sinfo->checks[i] != NULL) {
	    flag = sinfo->flags[i];
	    finfo->save_flags &= ~flag;
	    if (button_is_active(sinfo->checks[i])) {
		fname = entry_box_get_trimmed_text(sinfo->entries[i]);
		if (fname != NULL) {
		    finfo->save_flags |= flag;
		    reset_finfo_filename(finfo, i, fname);
		}
	    }
	}
    }

    sinfo->retval = 0;
    gtk_widget_destroy(sinfo->dialog);
}

static gchar *get_pkg_text_filename (function_info *finfo,
				     const char *pkgname,
				     const char **ids,
				     int i)
{
    const char *s;
    gchar *fname = NULL;

    s = function_package_get_string(finfo->pkg, ids[i]);

    if (s != NULL) {
	fname = g_strdup(s);
    } else if (i == 0) {
	fname = g_strdup_printf("%s_sample.inp", pkgname);
    } else if (i == 1) {
	fname = g_strdup_printf("%s_help.txt", pkgname);
    } else {
	fname = g_strdup_printf("%s_gui_help.txt", pkgname);
    }

    return fname;
}

static void sensitize_auxname_entry (GtkToggleButton *button,
				     GtkWidget *w)
{
    gboolean s = gtk_toggle_button_get_active(button);
	
    gtk_widget_set_sensitive(w, s);
}

static void nullify_spec_dialog (GtkWidget *w, function_info *finfo)
{
    finfo->specdlg = NULL;
}

static int gfn_spec_save_dialog (function_info *finfo,
				 const char **texts)
{
    const gchar *msgs[] = {
	N_("Save sample script as"),
	N_("Save help text as"),
	N_("Save GUI help as")
    };
    const char *ids[] = {
	"sample-fname",
	"help-fname",
	"gui-help-fname"
    };
    int flags[] = {
	WRITE_SAMPFILE,
	WRITE_HELPFILE,
	WRITE_GUI_HELP
    };    
    struct spec_info sinfo;
    GtkWidget *dialog, *entry;
    GtkWidget *vbox, *hbox, *w;
    GtkWidget *table;
    const char *pkgname;
    gchar *tmp;
    int i, j, n = 0;

    sinfo.retval = GRETL_CANCEL;
    sinfo.finfo = finfo;
    sinfo.flags = flags;

    finfo->specdlg = sinfo.dialog = dialog =
	gretl_dialog_new(NULL, finfo->dlg, GRETL_DLG_BLOCK);
    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(nullify_spec_dialog), finfo);
    
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    pkgname = function_package_get_name(finfo->pkg);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    tmp = g_strdup_printf(_("Saving %s.spec: also save ancillary file(s)?"),
			  pkgname);
    w = gtk_label_new(tmp);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);

    for (i=0; i<3; i++) {
	n += (texts[i] != NULL);
    }

    table = gtk_table_new(n, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);
    gtk_box_pack_start(GTK_BOX(vbox), table, FALSE, FALSE, 5);

    j = 0;
    for (i=0; i<3; i++) {
	if (texts[i] == NULL) {
	    sinfo.checks[i] = sinfo.entries[i] = NULL;
	    continue;
	}
	tmp = get_pkg_text_filename(finfo, pkgname, ids, i);
	w = sinfo.checks[i] = gtk_check_button_new_with_label(_(msgs[i]));
	if (finfo->save_flags & flags[i]) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
	}
	gtk_table_attach_defaults(GTK_TABLE(table), w, 0, 1, j, j+1);
	entry = sinfo.entries[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(entry), 64);
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 28);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	g_signal_connect(G_OBJECT(w), "toggled", 
			 G_CALLBACK(sensitize_auxname_entry),
			 entry);
	gtk_table_attach_defaults(GTK_TABLE(table), entry, 1, 2, j, j+1);
	g_free(tmp);
	j++;
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    cancel_delete_button(hbox, dialog);
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(spec_save_ok), &sinfo);
    gtk_widget_grab_default(w);

    gtk_widget_show_all(dialog);

    return sinfo.retval;
}

/* callback from file selector on saving package spec
   file: the default location should match that of the
   gfn file
*/

void get_gfn_dir (char *dirname, gpointer p)
{
    function_info *finfo = (function_info *) p;
    char *s;

    strcpy(dirname, finfo->fname);
    s = strrchr(dirname, SLASH);
    if (s != NULL) {
	*s = '\0';
    } else {
	*dirname = '\0';
    }
}

static void gfn_to_spec_callback (function_info *finfo)
{
    const char *texts[] = {
	finfo->sample,
	finfo->help,
	finfo->gui_help
    };
    int resp = 0;

    if (finfo->specdlg != NULL) {
	gtk_window_present(GTK_WINDOW(finfo->specdlg));
	return;
    }    
    
    if (finfo->pkg == NULL) {
	warnbox(_("Please save your package first"));
	return;
    }

    if (texts[0] == NULL) {
	texts[0] = function_package_get_string(finfo->pkg, "sample-script");
    }
    if (texts[1] != NULL && !strncmp(texts[1], "pdfdoc", 6)) {
	texts[1] = NULL;
    }    
    if (texts[2] == NULL) {
	texts[2] = function_package_get_string(finfo->pkg, "gui-help");
    }

    if (texts[0] != NULL || texts[1] != NULL || texts[2] != NULL) {
	resp = gfn_spec_save_dialog(finfo, texts);
    }

    if (!canceled(resp)) {
	file_selector_with_parent(SAVE_GFN_SPEC, FSEL_DATA_MISC, 
				  finfo, finfo->dlg);
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
	strings_array_free(names, n);
	names = NULL;
    }

    return names;
}

/* Convert from lists of functions by index numbers, @publist and
   @privlist, to arrays of public and private functions by name.
   Record in the @changed location whether or not there's any change
   relative to the existing lists of function names in @finfo.
*/

static int finfo_reset_function_names (function_info *finfo,
				       const int *pub_ids,
				       const int *priv_ids,
				       int *changed)
{
    char **pubnames = NULL;
    char **privnames = NULL;
    int npub = pub_ids[0];
    int npriv = (priv_ids == NULL)? 0 : priv_ids[0];
    int err = 0;

    *changed = 0;

    pubnames = get_function_names(pub_ids, &err);
    if (!err && npriv > 0) {
	privnames = get_function_names(priv_ids, &err);
    }

    if (!err) {
	if (npub != finfo->n_pub || npriv != finfo->n_priv) {
	    /* we know that something has changed */
	    *changed = 1;
	} else {
	    /* we'll have to check the arrays for any changes */
	    *changed = strings_array_cmp(pubnames, finfo->pubnames, npub);
	    if (*changed == 0 && npriv > 0) {
		*changed = strings_array_cmp(privnames, finfo->privnames,
					     npriv);
	    }
	}
    }

    if (*changed == 0) {
	/* trash the new function-name arrays */
	strings_array_free(pubnames, npub);
	strings_array_free(privnames, npriv);
    } else if (!err) {
	/* replace the old function-name arrays */
	strings_array_free(finfo->pubnames, finfo->n_pub);
	finfo->pubnames = pubnames;
	finfo->n_pub = npub;
	strings_array_free(finfo->privnames, finfo->n_priv);
	finfo->privnames = privnames;
	finfo->n_priv = npriv;
	finfo->active = finfo->pubnames[0];
    }

    return err;
}

static int finfo_set_function_names (function_info *finfo,
				     const int *publist,
				     const int *privlist)
{
    int npriv = (privlist == NULL)? 0 : privlist[0];
    int err = 0;

    finfo->pubnames = get_function_names(publist, &err);
    if (!err) {
	finfo->n_pub = publist[0];
    }

    if (!err && npriv > 0) {
	finfo->privnames = get_function_names(privlist, &err);
	if (!err) {
	    finfo->n_priv = npriv;
	}
    }

    return err;
}

static void func_selector_set_strings (function_info *finfo, 
				       GtkWidget *ifmenu)
{
    int i, n = 0;

    for (i=0; i<finfo->n_pub; i++) {
	combo_box_append_text(ifmenu, finfo->pubnames[i]);
	n++;
    }

    for (i=0; i<finfo->n_priv; i++) {
	gchar *s = g_strdup_printf("%s (%s)", finfo->privnames[i], _("private"));

	combo_box_append_text(ifmenu, s);
	g_free(s);
	n++;
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(ifmenu), 0);
    gtk_widget_set_sensitive(ifmenu, n > 1);
}

static GtkWidget *active_func_selector (function_info *finfo)
{
    GtkWidget *ifmenu = gtk_combo_box_text_new();

    func_selector_set_strings(finfo, ifmenu);
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
	N_("Panel data"),
	N_("No dataset needed")
    };
    GtkWidget *datamenu, *tmp;
    int j;

    tmp = gtk_label_new(_("Data requirement"));
    gtk_misc_set_alignment(GTK_MISC(tmp), 1.0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
    gtk_widget_show(tmp);

    datamenu = gtk_combo_box_text_new();
    for (j=0; j<=FN_NODATA_OK; j++) {
	combo_box_append_text(datamenu, _(datareq[j]));
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
    *pl = v - *maj * 10000 - *min * 100;
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
    int x, y, z;
    int ymin = 0;

    get_maj_min_pl(finfo->minver, &maj, &min, &pl);
    sscanf(GRETL_VERSION, "%d.%d.%d", &x, &y, &z);
    if (x == 1) {
	ymin = 8;
    }

    tmp = gtk_label_new(_("Minimum gretl version"));
    gtk_misc_set_alignment(GTK_MISC(tmp), 1.0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
    gtk_widget_show(tmp);

    hbox = gtk_hbox_new(FALSE, 0);

    /* major version number */
    spin = gtk_spin_button_new_with_range(1, x, 1);
    if (maj > 1) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), (double) maj);
    }
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 2);
    g_object_set_data(G_OBJECT(spin), "level", GINT_TO_POINTER(1));
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(adjust_minver), finfo);
    tmp = gtk_label_new(".");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 2);

    /* minor version number */
    spin = gtk_spin_button_new_with_range(ymin, y, 1);
    if (min > 0) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), (double) min);
    }
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 2);
    g_object_set_data(G_OBJECT(spin), "level", GINT_TO_POINTER(2));
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(adjust_minver), finfo);
    tmp = gtk_label_new(".");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 2);

    /* "patch-level" */
    spin = gtk_spin_button_new_with_range(0, z + 1, 1);
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

    return TRUE;
}

static gint query_save_package (GtkWidget *w, GdkEvent *event, 
				function_info *finfo)
{
    if (finfo->modified) {
	int resp =
	    yes_no_cancel_dialog("gretl", _("Save changes?"), w);

	if (resp == GRETL_CANCEL) {
	    return TRUE;
	} else if (resp == GRETL_YES) {
	    return finfo_save(finfo);
	}
    }

    return FALSE;
}

static void check_pkg_callback (GtkWidget *widget, function_info *finfo)
{
    validate_package_file(finfo->fname, 1);
}

static GtkTreeStore *make_menu_attachment_tree (function_info *finfo,
						GtkTreePath **ppath,
						int modelwin)
{
    const char *main_items = 
	"0 Tools\n"
	"0 Data\n"
	"0 View\n"
	"1 GraphVars\n"
	"2 MultiPlots\n"
	"0 Add\n"
	"0 Sample\n"
	"0 Variable\n"
	"1 URTests\n"
	"1 Filter\n"
	"0 Model\n"
	"1 ivreg\n"  
	"1 LinearModels\n"
	"1 LimdepModels\n"
	"2 logit\n"
	"2 probit\n"
	"1 TSModels\n"
	"2 CointMenu\n"
	"1 PanelModels\n"
	"1 RobustModels\n";
    const char *model_items = 
	"0 Tests\n"
	"0 Save\n"
	"0 Graphs\n"
	"0 Analysis\n";
    const char *leaders[] = {
	"MAINWIN",
	"MODELWIN"
    };
    const char *leader;
    GtkTreeStore *store;
    GtkTreeIter iter;
    GtkTreeIter parents[2];
    GtkTreeIter *iterp;
    gchar *path, *ustr;
    const char *s;
    char words[3][16];
    char *word;
    int level;

    store = gtk_tree_store_new(2, G_TYPE_STRING, G_TYPE_STRING);

    leader = modelwin ? leaders[1] : leaders[0];
    s = modelwin ? model_items : main_items;

    while (*s) {
	level = atoi(s);
	word = words[level];
	s += 2;
	sscanf(s, "%s", word);
	ustr = user_friendly_menu_path(word, modelwin);
	if (level == 0) {
	    path = g_strdup_printf("%s/%s", leader, words[0]);
	    gtk_tree_store_append(store, &parents[0], NULL);
	    iterp = &parents[0];
	} else if (level == 1) {
	    path = g_strdup_printf("%s/%s/%s", leader, words[0], 
				   words[1]);
	    gtk_tree_store_append(store, &parents[1], &parents[0]);
	    iterp = &parents[1];
	} else {
	    path = g_strdup_printf("%s/%s/%s/%s", leader, words[0], 
				   words[1], words[2]); 
	    gtk_tree_store_append(store, &iter, &parents[1]);
	    iterp = &iter;
	}
	gtk_tree_store_set(store, iterp, 0, ustr, 1, path, -1);
	if (finfo->menupath != NULL && !strcmp(path, finfo->menupath)) {
	    /* record the path of pre-selected menu entry */
	    *ppath = gtk_tree_model_get_path(GTK_TREE_MODEL(store), iterp);
	}	
	g_free(path);
	g_free(ustr);
	s += strlen(word) + 1;
    }

    return store;
}

static GtkWidget *add_menu_navigator (GtkWidget *holder, 
				      function_info *finfo,
				      int modelwin)
{
    GtkTreeStore *store;
    GtkWidget *view;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    GtkTreePath *path = NULL;

    store = make_menu_attachment_tree(finfo, &path, modelwin);
    if (store == NULL) {
	return NULL;
    }

    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);
    g_object_set(view, "enable-tree-lines", TRUE, NULL);
    g_object_ref(G_OBJECT(view));

    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes("",
						      renderer,
						      "text", 0,
						      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);

    if (path != NULL) {
	gtk_tree_view_expand_to_path(GTK_TREE_VIEW(view), path);
	gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(view), path,
				     NULL, FALSE, 0, 0);
	gtk_tree_selection_select_path(select, path);
	gtk_tree_path_free(path);
    }

    if (finfo->treewin == NULL) {
	GtkWidget *sw;

	sw = finfo->treewin = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				       GTK_POLICY_AUTOMATIC, 
				       GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					    GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(holder), sw);
	gtk_widget_set_size_request(sw, 150, -1);
    }

    if ((modelwin == 0 && finfo->menuwin != MODEL_WINDOW) ||
	(modelwin == 1 && finfo->menuwin == MODEL_WINDOW)) {
	gtk_container_add(GTK_CONTAINER(finfo->treewin), view);
	finfo->currtree = view;
    } else {
	finfo->alttree = view;
    }

    gtk_tree_view_columns_autosize(GTK_TREE_VIEW(view));

    return view;
}

static void switch_menu_view (GtkComboBox *combo,
			      function_info *finfo) 
{
    int w = gtk_combo_box_get_active(combo);
    GtkWidget *sw = finfo->treewin;

    if (w == NO_WINDOW) {
	if (finfo->currtree != NULL) {
	    gtk_widget_set_sensitive(finfo->currtree, FALSE);
	}
    } else if (w == MAIN_WINDOW) {
	if (finfo->currtree == finfo->modeltree) {
	    gtk_container_remove(GTK_CONTAINER(sw), finfo->modeltree);
	    finfo->alttree = finfo->modeltree;
	    gtk_widget_show(finfo->maintree);
	    gtk_container_add(GTK_CONTAINER(sw), finfo->maintree);
	    finfo->currtree = finfo->maintree;
	}
	gtk_widget_set_sensitive(finfo->currtree, TRUE);
    } else if (w == MODEL_WINDOW) {
	if (finfo->currtree == finfo->maintree) {
	    gtk_container_remove(GTK_CONTAINER(sw), finfo->maintree);
	    finfo->alttree = finfo->maintree;
	    gtk_widget_show(finfo->modeltree);
	    gtk_container_add(GTK_CONTAINER(sw), finfo->modeltree);
	    finfo->currtree = finfo->modeltree;
	}
	gtk_widget_set_sensitive(finfo->currtree, TRUE);
    }	
}

static void add_data_files_entries (GtkWidget *holder,
				    function_info *finfo)
{
    const char *msg = N_("You may add or delete names of data"
			 "files to be included in the package.");
    GtkWidget *w, *hbox, *entry;
    int i;

    w = gtk_label_new(_(msg));
    gtk_label_set_line_wrap(GTK_LABEL(w), TRUE);
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(holder), hbox, FALSE, FALSE, 5);

    for (i=0; i<N_FILE_ENTRIES; i++) {
	hbox = gtk_hbox_new(FALSE, 5);
	finfo->file_entries[i] = entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 32);
	if (i < finfo->n_files) {
	    gtk_entry_set_text(GTK_ENTRY(entry), finfo->datafiles[i]);
	}
	gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(holder), hbox, FALSE, FALSE, 5);
    }
}

static void gui_help_text_callback (GtkButton *b, function_info *finfo)
{
    const char *pkgname;
    gchar *title;
    PRN *prn = NULL;

    if (finfo->gui_helpwin != NULL) {
	gtk_window_present(GTK_WINDOW(finfo->gui_helpwin->main));
	return;
    }

    if (finfo->gui_help == NULL) {
	const char *msg = 
	    N_("This package has no GUI-specific help text at present.\n"
	       "Would you like to add some?");

	if (yes_no_dialog(NULL, _(msg), finfo->extra) != GRETL_YES) {
	    return;
	}
    }

    if (bufopen(&prn)) {
	return;
    }

    pkgname = function_package_get_name(finfo->pkg);
    title = g_strdup_printf("%s gui-help", pkgname);

    if (finfo->gui_help != NULL) {
	pputs(prn, finfo->gui_help);
	pputc(prn, '\n');
    } 

    finfo->gui_helpwin = view_buffer(prn, 72, 350, title,
				     EDIT_PKG_GHLP, finfo);
    g_object_set_data(G_OBJECT(finfo->gui_helpwin->main), "finfo",
		      finfo);
    g_signal_connect(G_OBJECT(finfo->gui_helpwin->main), "destroy",
		     G_CALLBACK(nullify_gui_helpwin), finfo);
    g_free(title);
}

static void add_menu_attach_top (GtkWidget *holder,
				 function_info *finfo)
{
    GtkWidget *w, *hbox, *entry;
    GtkWidget *combo;

    /* menu label entry */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Label"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), 36);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 32);
    if (finfo->menulabel != NULL) {
	gtk_entry_set_text(GTK_ENTRY(entry), finfo->menulabel);
    }
    gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 5);
    g_object_set_data(G_OBJECT(finfo->extra), "label-entry", entry);
    gtk_box_pack_start(GTK_BOX(holder), hbox, FALSE, FALSE, 5);

    /* menu attachment combo */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Window"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    combo = gtk_combo_box_text_new();
    combo_box_append_text(combo, _("none"));    
    combo_box_append_text(combo, _("main window"));    
    combo_box_append_text(combo, _("model window"));    
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), finfo->menuwin);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(switch_menu_view), finfo);

    /* gui-help button */
    w = gtk_button_new_with_label("GUI help text");
    g_signal_connect(G_OBJECT(w), "clicked",
		     G_CALLBACK(gui_help_text_callback), finfo);
    gtk_box_pack_end(GTK_BOX(hbox), w, FALSE, FALSE, 5);

    /* complete the packing */
    gtk_box_pack_start(GTK_BOX(holder), hbox, FALSE, FALSE, 5);
}

/* pertaining to the "extra properties" dialog: check for
   any changes in relation to menu attachment
*/

static int process_menu_attachment (function_info *finfo,
				    gboolean make_changes,
				    int *focus_label)
{
    GtkTreeSelection *selection;
    GtkTreeModel *model;
    GtkTreeIter iter;
    GtkWidget *view, *entry;
    gchar *label;
    int changed = 0;

    view = finfo->currtree;
    
    entry = g_object_get_data(G_OBJECT(finfo->extra), "label-entry");
    label = entry_box_get_trimmed_text(entry);

    if (label == NULL || *label == '\0') {
	if (finfo->menulabel != NULL) {
	    if (make_changes) {
		g_free(finfo->menulabel);
		finfo->menulabel = NULL;
	    }
	    changed = 1;
	}
    } else if (finfo->menulabel == NULL || strcmp(finfo->menulabel, label)) {
	if (make_changes) {
	    g_free(finfo->menulabel);
	    finfo->menulabel = label;
	    label = NULL;
	}
	changed = 1;
    }

    g_free(label);

    if (!gtk_widget_is_sensitive(view)) {
	/* no menu attachment at present */
	if (finfo->menupath != NULL) {
	    if (make_changes) {
		g_free(finfo->menupath);
		finfo->menupath = NULL;
	    }
	    changed = 1;
	}
	finfo->menuwin = NO_WINDOW;
	goto finish;
    }	

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));

    if (!gtk_tree_selection_get_selected(selection, &model, &iter)) {
	if (finfo->menupath != NULL) {
	    if (make_changes) {
		g_free(finfo->menupath);
		finfo->menupath = NULL;
	    }
	    changed = 1;
	}
    } else {
	gchar *newpath = NULL;

	gtk_tree_model_get(model, &iter, 1, &newpath, -1);
	if (finfo->menupath == NULL || strcmp(finfo->menupath, newpath)) {
	    if (make_changes) {
		g_free(finfo->menupath);
		finfo->menupath = newpath;
		newpath = NULL;
	    }
	    changed = 1;
	}
	g_free(newpath);
    }

    finfo_set_menuwin(finfo);

 finish:

    if (make_changes) {
	if (finfo->menupath != NULL && *finfo->menupath != '\0' &&
	    (finfo->menulabel == NULL || *finfo->menulabel == '\0')) {
	    warnbox(_("To create a menu attachment, you must supply a label."));
	    *focus_label = 1;
	}
    }

    return changed;
}

static int want_no_print_toggle (int role)
{
    return role != UFUN_GUI_PRECHECK &&
	role != UFUN_BUNDLE_PRINT;
}

/* pertaining to the "extra properties" dialog: check for
   any changes in relation to the special functions table
*/

static int process_special_functions (function_info *finfo,
				      gboolean make_changes)
{
    GtkWidget **c_array;
    const char *oldfun;
    gchar *newfun;
    int n_changed = 0;
    int i, err = 0;

    c_array = g_object_get_data(G_OBJECT(finfo->extra), "combo-array");

    /* For each special function slot, check to see if the 
       currently selected function differs from what was
       present originally, and if so update the record in
       @finfo->specials. The changes are not yet saved to
       the function package itself.
    */

    for (i=0; i<N_SPECIALS && !err; i++) {
	int role = i + 1;
	
	if (gtk_widget_is_sensitive(c_array[i])) {
	    int fn_changed = 0, attr_changed = 0;
	    int newnull = 0, oldnull = 0;
	    unsigned char attr = 0;
	    GtkWidget *cb;

	    /* retrieve and check the selected name */
	    newfun = combo_box_get_active_text(GTK_COMBO_BOX(c_array[i]));
	    newnull = (newfun == NULL || *newfun == '\0' || 
		       !strcmp(newfun, "none"));

	    /* retrieve and check what was there before */
	    oldfun = finfo->specials[i];
	    oldnull = (oldfun == NULL || *oldfun == '\0' || 
		       !strcmp(oldfun, "none"));

	    if (want_no_print_toggle(role)) {
		/* retrieve the no-print attribute? */
		cb = g_object_get_data(G_OBJECT(c_array[i]), "np-toggle");
		if (cb != NULL && button_is_active(cb)) {
		    attr |= UFUN_NOPRINT;
		}
	    }

	    if (role == UFUN_GUI_MAIN) {
		/* retrieve the menu-only attribute? */
		cb = g_object_get_data(G_OBJECT(c_array[i]), "mo-toggle");
		if (cb != NULL && button_is_active(cb)) {
		    attr |= UFUN_MENU_ONLY;
		}
	    }

	    if (oldnull && !newnull) {
		fn_changed = 1;
	    } else if (!oldnull && newnull) {
		fn_changed = 1;
	    } else if (!oldnull && !newnull) {
		fn_changed = strcmp(newfun, oldfun);
	    }

	    if (attr != finfo->gui_attrs[i]) {
		attr_changed = 1;
	    }

	    if (make_changes) {
		if (fn_changed) {
		    free(finfo->specials[i]);
		    finfo->specials[i] = gretl_strdup(newfun);
		}
		if (attr_changed) {
		    finfo->gui_attrs[i] = attr;
		}
	    }

	    if (fn_changed || attr_changed) {
		n_changed++;
	    }
	    
	    g_free(newfun);
	}
    }

    return n_changed;
}

#define must_be_private(r) (r == UFUN_GUI_PRECHECK)
#define must_be_public(r) (r != UFUN_GUI_PRECHECK)

/* After adding or deleting functions, check that any
   selected "specials" are still valid: the selected
   funtion has not been removed from the package, nor
   has its public/private status been changed such as
   to disqualify it from playing the given role. If a
   selection has been invalidated, null it out.
*/

static void verify_selected_specials (function_info *finfo)
{
    const char *seek;
    int i, j, found, role;

    for (i=0; i<N_SPECIALS; i++) {
	role = i + 1;
	if (finfo->specials[i] != NULL) {
	    /* a selection was made */
	    seek = finfo->specials[i];
	    found = 0;
	    if (!must_be_private(role)) {
		/* try the public interface list */
		for (j=0; j<finfo->n_pub && !found; j++) {
		    if (!strcmp(seek, finfo->pubnames[j])) {
			found = 1;
		    }
		}
	    }
	    if (!found && !must_be_public(role)) {
		/* try the private interface list */
		for (j=0; j<finfo->n_priv && !found; j++) {
		    if (!strcmp(seek, finfo->privnames[j])) {
			found = 1;
		    }
		}
	    }
	    if (!found) {
		/* gone bad */
		free(finfo->specials[i]);
		finfo->specials[i] = NULL;
	    }
	}
    }
}

static int data_file_check_existence (function_info *finfo,
				      const char *fname)
{
    char *p, test[FILENAME_MAX];
    struct stat buf = {0};

    strcpy(test, finfo->fname);
    p = strrchr(test, SLASH);
    if (p != NULL) {
	*p = '\0';
	strcat(p, fname);
    } else {
	strcpy(test, fname);
    }

    if (gretl_stat(test, &buf) != 0) {
	gchar *msg;

	msg = g_strdup_printf(_("Couldn't find %s"), test);
	msgbox(msg, GTK_MESSAGE_WARNING, finfo->extra);	
	g_free(msg);
	return 1;
    }

    return 0;
}

/* pertaining to the "extra properties" dialog: check for
   any changes in relation to included data files
*/

static int process_data_file_names (function_info *finfo,
				    gboolean make_changes)
{
    gchar *fname;
    int i, nf = 0;
    int changed = 0;

    for (i=0; i<N_FILE_ENTRIES; i++) {
	fname = entry_box_get_trimmed_text(finfo->file_entries[i]);
	if (fname != NULL) {
	    nf++;
	    if (i < finfo->n_files &&
		strcmp(fname, finfo->datafiles[i])) {
		changed = 1;
	    }
	}
	g_free(fname);
    }

    if (!changed && nf != finfo->n_files) {
	/* added or deleted */
	changed = 1;
    }

    if (changed && make_changes) {
	strings_array_free(finfo->datafiles, finfo->n_files);
	if (nf == 0) {
	    finfo->datafiles = NULL;
	    finfo->n_files = 0;
	} else {
	    finfo->datafiles = strings_array_new(nf);
	    if (finfo->datafiles != NULL) {
		finfo->n_files = nf;
	    }
	}

	if (finfo->datafiles != NULL) {
	    int j = 0, err = 0;
	    
	    for (i=0; i<N_FILE_ENTRIES; i++) {
		fname = entry_box_get_trimmed_text(finfo->file_entries[i]);
		if (fname != NULL) {
		    if (make_changes && err == 0) {
			err = data_file_check_existence(finfo, fname);
		    }
		    finfo->datafiles[j++] = gretl_strdup(fname);
		}
		g_free(fname);
	    }
	}	
    }

    return changed;
}

static int process_extra_properties (function_info *finfo,
				     gboolean make_changes)
{
    int focus_label = 0;
    int changed = 0;

    changed += process_special_functions(finfo, make_changes);

    changed += process_menu_attachment(finfo, make_changes, &focus_label);
    if (focus_label) {
	GtkWidget *w = g_object_get_data(G_OBJECT(finfo->extra), "label-entry");

	gtk_widget_grab_focus(w);
    }

    changed += process_data_file_names(finfo, make_changes);

    if (changed && make_changes) {
	finfo_set_modified(finfo, TRUE);
    }

    return changed;
}

static void extra_properties_apply (GtkWidget *w, function_info *finfo)
{
    process_extra_properties(finfo, TRUE);
}

static void extra_properties_close (GtkWidget *w, function_info *finfo)
{
    int changed = process_extra_properties(finfo, FALSE);

    if (changed) {
	int resp = yes_no_cancel_dialog(NULL, _("Apply changes?"),
					finfo->extra);

	if (resp == GRETL_CANCEL) {
	    return;
	} else if (resp == GRETL_YES) {
	    process_extra_properties(finfo, TRUE);
	}
    }
    
    gtk_widget_destroy(finfo->extra);
}

static gint query_save_extra_props (GtkWidget *w, GdkEvent *event, 
				    function_info *finfo)
{
    int changed = process_extra_properties(finfo, FALSE);

    if (changed) {
	int resp = yes_no_cancel_dialog(NULL, _("Apply changes?"),
					finfo->extra);

	if (resp == GRETL_CANCEL) {
	    return TRUE;
	} else if (resp == GRETL_YES) {
	    process_extra_properties(finfo, TRUE);
	}
    }
    
    return FALSE;
}

static void sensitize_attr_toggles (GObject *obj, gboolean s)
{
    GtkWidget *cb;

    cb = g_object_get_data(obj, "np-toggle");
    if (cb != NULL) {
	gtk_widget_set_sensitive(cb, s);
    }

    cb = g_object_get_data(obj, "mo-toggle");
    if (cb != NULL) {
	gtk_widget_set_sensitive(cb, s);
    }    
}

/* Prevent the user from assigning a given function to more
   then one special role: when a selection is changed, if
   the given function is already selected for a different
   role, deselect it in that role.
*/

static void special_changed_callback (GtkComboBox *this, 
				      function_info *finfo)
{
    GtkWidget **c_array;
    GtkWidget *other;
    gchar *s0, *si;
    int i, dup = 0;

    if (gtk_combo_box_get_active(this) == 0) {
	/* selected "none" */
	sensitize_attr_toggles(G_OBJECT(this), FALSE);
	return;
    }

    sensitize_attr_toggles(G_OBJECT(this), TRUE);
    s0 = combo_box_get_active_text(this);
    c_array = g_object_get_data(G_OBJECT(finfo->extra), "combo-array");

    for (i=0; i<N_SPECIALS && !dup; i++) {
	other = c_array[i];
	if (other != GTK_WIDGET(this)) {
	    si = combo_box_get_active_text(GTK_COMBO_BOX(other));
	    if (!strcmp(si, s0)) {
		/* switch to "none" */
		gtk_combo_box_set_active(GTK_COMBO_BOX(other), 0);
		sensitize_attr_toggles(G_OBJECT(other), FALSE);
		dup = 1;
	    }
	    g_free(si);
	}
    }

    g_free(s0);
}

static void finfo_extra_help (GtkWidget *w, function_info *finfo)
{
    GtkWidget *notebook;
    gint page;

    notebook = g_object_get_data(G_OBJECT(finfo->extra), "book");
    page = gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));
    
    if (page == 0) {
	show_gui_help(GUI_FUNCS);
    } else if (page == 1) {
	show_gui_help(MENU_ATTACH);
    } else {
	show_gui_help(PKG_FILES);
    }
}

static void unref_trees (GtkWidget *w, function_info *finfo)
{
    g_object_unref(G_OBJECT(finfo->maintree));
    g_object_unref(G_OBJECT(finfo->modeltree));
}

/* The following function supports three "notebook" tabs. The first
   allows the user to select functions in the package for the various
   "special" package roles (e.g. gui-main, bundle-print). The second
   allows for selection of a menu attachment point and GUI label. The
   third allows for specification of additional data to be included in
   the package.
*/

static void extra_properties_dialog (GtkWidget *w, function_info *finfo)
{
    GtkWidget *dlg, *combo, *table;
    GtkWidget *tmp, *vbox, *hbox;
    GtkWidget **combo_array;
    GtkWidget *notebook;
    const char *funname;
    const char *key;
    const char *special;
    int tabcols = 4;
    int nfuns, i, j;

    if (finfo->extra != NULL) {
	gtk_window_present(GTK_WINDOW(finfo->extra));
	return;
    }

    if (finfo->pkg == NULL) {
	warnbox(_("Please save your package first"));
	return;
    }

    finfo->maintree = finfo->modeltree = NULL;
    finfo->currtree = finfo->alttree = NULL;
    finfo->treewin = NULL;

    dlg = gretl_dialog_new(_("gretl: extra properties"), finfo->dlg,
			   GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);
    finfo->extra = dlg;
    g_signal_connect(G_OBJECT(dlg), "delete-event",
		     G_CALLBACK(query_save_extra_props), finfo);    
    g_signal_connect(G_OBJECT(dlg), "destroy",
		     G_CALLBACK(gtk_widget_destroyed), &finfo->extra);
    g_signal_connect(G_OBJECT(dlg), "destroy",
		     G_CALLBACK(unref_trees), finfo);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(vbox), notebook, TRUE, TRUE, 0);
    g_object_set_data(G_OBJECT(dlg), "book", notebook);

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);
    tmp = gtk_label_new(_("Special functions"));
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, tmp);  

    table = gtk_table_new(N_SPECIALS, tabcols, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);

    nfuns = finfo->n_priv + finfo->n_pub;
    combo_array = g_malloc(N_SPECIALS * sizeof *combo_array);
    g_object_set_data_full(G_OBJECT(dlg), "combo-array", 
			   combo_array, g_free);

    /* For each "special" function role, test the functions
       in finfo->pkg to see if they qualify as candidates for
       that role; if so, add them to the combo selector.
    */

    for (i=0; i<N_SPECIALS; i++) {
	int n_cands = 0;
	int selected = 0;
	int role = i + 1;

	key = package_role_get_key(role);
	special = finfo->specials[i];
	combo = gtk_combo_box_text_new();
	combo_box_append_text(combo, "none");

	for (j=0; j<nfuns; j++) {
	    if (j < finfo->n_priv && !must_be_public(role)) {
		funname = finfo->privnames[j];
	    } else if (j >= finfo->n_priv && !must_be_private(role)) {
		funname = finfo->pubnames[j - finfo->n_priv];
	    } else {
		continue;
	    }
	    if (function_ok_for_package_role(funname, role)) {
		combo_box_append_text(combo, funname);
		if (special != NULL && !selected && !strcmp(special, funname)) {
		    selected = n_cands + 1;
		}
		n_cands++;
	    }
	}

	tmp = gtk_label_new(key);
	gtk_misc_set_alignment(GTK_MISC(tmp), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(table), tmp, 0, 1, i, i+1);
	gtk_table_attach_defaults(GTK_TABLE(table), combo, 1, 2, i, i+1);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), selected);
	if (n_cands == 0) {
	    gtk_widget_set_sensitive(combo, FALSE);
	} else {
	    g_signal_connect(G_OBJECT(combo), "changed", 
			     G_CALLBACK(special_changed_callback), 
			     finfo);
	}
	combo_array[i] = combo;

	if (want_no_print_toggle(role)) {
	    GtkWidget *cb = gtk_check_button_new_with_label("no-print");

	    gtk_table_attach_defaults(GTK_TABLE(table), cb, 2, 3, i, i+1);
	    g_object_set_data(G_OBJECT(combo), "np-toggle", cb);
	    gtk_widget_set_sensitive(cb, selected > 0);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(cb),
					 finfo->gui_attrs[i] & UFUN_NOPRINT);
	}

	if (role == UFUN_GUI_MAIN) {
	    GtkWidget *cb = gtk_check_button_new_with_label("menu-only");

	    gtk_table_attach_defaults(GTK_TABLE(table), cb, 3, 4, i, i+1);
	    g_object_set_data(G_OBJECT(combo), "mo-toggle", cb);
	    gtk_widget_set_sensitive(cb, selected > 0);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(cb),
					 finfo->gui_attrs[i] & UFUN_MENU_ONLY);
	}	
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), table, FALSE, FALSE, 5);

    /* the menu attachment page */

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);
    tmp = gtk_label_new(_("Menu attachment"));
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, tmp);

    add_menu_attach_top(vbox, finfo);
    finfo->maintree = add_menu_navigator(vbox, finfo, 0);
    finfo->modeltree = add_menu_navigator(vbox, finfo, 1);

    if (finfo->menuwin == NO_WINDOW) {
	gtk_widget_set_sensitive(finfo->currtree, FALSE);
    }

    /* the data files page */
    
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);
    tmp = gtk_label_new(_("Data files"));
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, tmp);
    add_data_files_entries(vbox, finfo);

    /* the common buttons area */
    
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Apply button */
    tmp = apply_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(extra_properties_apply), finfo);
    gtk_widget_grab_default(tmp);

    /* Close button */
    tmp = close_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
                     G_CALLBACK(extra_properties_close), finfo);

    /* Help button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_HELP);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       tmp, TRUE);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(finfo_extra_help),
		     finfo);

    gtk_widget_show_all(dlg);
}

int package_editor_exit_check (GtkWidget *w)
{
    function_info *finfo;

    finfo = g_object_get_data(G_OBJECT(w), "finfo");
    
    if (finfo != NULL && finfo->modified) {
	gtk_window_present(GTK_WINDOW(w));
	return query_save_package(w, NULL, finfo);
    }

    return FALSE;
}

/* return non-zero if @w is the window of an editor working
   on @pkgname
*/

int query_package_editor (GtkWidget *w, const char *pkgname)
{
    function_info *finfo;

    finfo = g_object_get_data(G_OBJECT(w), "finfo");
    
    if (finfo != NULL && finfo->pkg != NULL) {
	const char *myname = function_package_get_name(finfo->pkg);
	
	return strcmp(pkgname, myname) == 0;
    }

    return FALSE;
}

static void delete_dlg_callback (GtkWidget *button, function_info *finfo) 
{
    gint resp = 0;

    if (finfo->modified) {
	resp = query_save_package(finfo->dlg, NULL, finfo);
    }

    if (!resp) {
	gtk_widget_destroy(finfo->dlg); 
    }
}

/* Dialog for editing function package.  The user can get here in
   either of two ways: after selecting functions to put into a
   newly created package, or upon selecting an existing package
   for editing.
*/

static void finfo_dialog (function_info *finfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *vbox, *hbox;
    GtkTextBuffer *hbuf = NULL;
    const char *entry_labels[] = {
	N_("Author"),
	N_("Email"),
	N_("Version"),
	N_("Date (YYYY-MM-DD)"),
	N_("Package description")
    };
    char *entry_texts[] = {
	finfo->author,
	finfo->email,
	finfo->version,
	finfo->date,
	finfo->pkgdesc
    };
    int focused = 0;
    int i;

    finfo->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_default_size(GTK_WINDOW(finfo->dlg), 640, 500);

    if (finfo->fname != NULL) {
	gchar *title = title_from_filename(finfo->fname, EDIT_PKG_CODE, TRUE);

	gtk_window_set_title(GTK_WINDOW(finfo->dlg), title);
	g_free(title);
    } else {
	gtk_window_set_title(GTK_WINDOW(finfo->dlg), 
			     _("gretl: function package editor"));
    } 

    g_object_set_data(G_OBJECT(finfo->dlg), "finfo", finfo);
    gtk_widget_set_name(finfo->dlg, "pkg-editor");
    g_signal_connect(G_OBJECT(finfo->dlg), "delete-event",
		     G_CALLBACK(query_save_package), finfo);
    g_signal_connect(G_OBJECT(finfo->dlg), "destroy", 
		     G_CALLBACK(finfo_destroy), finfo);

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 5);
    gtk_container_add(GTK_CONTAINER(finfo->dlg), vbox);

    tbl = gtk_table_new(NENTRIES + 1, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 4);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 5);

    for (i=0; i<NENTRIES; i++) {
	GtkWidget *entry;

	label = gtk_label_new(_(entry_labels[i]));
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, i, i+1);

	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 40);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);

	finfo->entries[i] = entry;

	if (entry_texts[i] != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(entry), entry_texts[i]);
	    if (i == 3) {
		g_signal_connect(G_OBJECT(entry), "button-press-event",
				 G_CALLBACK(today_popup), &finfo->popup);
	    }
	} else if (i == 0) {
	    /* author */
	    const gchar *s = get_user_string();

	    if (s != NULL) {
		gtk_entry_set_text(GTK_ENTRY(entry), s);
	    }
	} else if (i == 1) {
	    /* email */
	    gtk_entry_set_text(GTK_ENTRY(entry), get_author_mail());
	} else if (i == 2) {
	    /* version */
	    gtk_entry_set_text(GTK_ENTRY(entry), "1.0");
	} else if (i == 3) {
	    /* date */
	    gtk_entry_set_text(GTK_ENTRY(entry), print_today());
	}

	if (i == 0 && entry_texts[i] == NULL) {
	    /* no author's name */
	    gtk_widget_grab_focus(entry);
	    focused = 1;
	} else if (i == 1 && !focused &&
		   (entry_texts[i] == NULL || *entry_texts[i] == '\0')) {
	    /* no email address */
	    gtk_widget_grab_focus(entry);
	    focused = 1;
	} else if (i == 2 && !focused) {
	    /* version number */
	    gtk_widget_grab_focus(entry);
	}

	g_signal_connect(GTK_EDITABLE(entry), "changed",
			 G_CALLBACK(pkg_changed), finfo);
    }

    add_minver_selector(tbl, i++, finfo);
    add_data_requirement_menu(tbl, i, finfo);

    hbox = label_hbox(vbox, _("Help text:"));

    finfo->text = editable_text_box(&hbuf);
    text_table_setup(vbox, finfo->text);

    if (finfo->help != NULL) {
	textview_set_text(finfo->text, finfo->help);
    }
    g_signal_connect(G_OBJECT(hbuf), "changed", 
		     G_CALLBACK(pkg_changed), finfo);

    /* table for buttons arrayed at foot of window */
    hbox = gtk_hbox_new(FALSE, 0);
    tbl = gtk_table_new(2, 4, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 4);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 4);
    gtk_table_set_col_spacing(GTK_TABLE(tbl), 2, 32);
    gtk_box_pack_start(GTK_BOX(hbox), tbl, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 2);

    /* first button row */

    /* 1: edit code button */
    button = gtk_button_new_with_label(_("Edit function code"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 0, 1, 0, 1);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_code_callback), finfo);

    /* 2: interface selector */
    finfo->codesel = active_func_selector(finfo);
    gtk_table_attach_defaults(GTK_TABLE(tbl), finfo->codesel, 1, 2, 0, 1);
    g_signal_connect(G_OBJECT(finfo->codesel), "changed",
		     G_CALLBACK(update_active_func), finfo);

    update_active_func(NULL, finfo);

    /* 3: extra package properties button */
    button = gtk_button_new_with_label(_("Extra properties"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 2, 3, 0, 1);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(extra_properties_dialog), finfo);

    /* 4: save-menu button */
    button = gtk_button_new_with_label(_("Save..."));
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 3, 4, 0, 1);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(pkg_save_popup), finfo);    

    /* second button row */

    /* 1: edit sample script button */
    button = gtk_button_new_with_label(_("Edit sample script"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 0, 1, 1, 2);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_sample_callback), finfo);

    /* 2: add/remove functions button */
    button = gtk_button_new_with_label(_("Add/remove functions"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 1, 2, 1, 2);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(add_remove_callback), finfo);

    /* 3: validate button */
    button = gtk_button_new_with_label(_("Validate"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 2, 3, 1, 2);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(check_pkg_callback), finfo);
    gtk_widget_set_sensitive(button, finfo->fname != NULL);
    finfo->validate = button;

    /* 4: close button */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 3, 4, 1, 2);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_dlg_callback), finfo);       

    finfo_set_modified(finfo, finfo->fname == NULL);

    window_list_add(finfo->dlg, SAVE_FUNCTIONS);
    gtk_widget_show_all(finfo->dlg);
}

static void web_get_login (GtkWidget *w, gpointer p)
{
    browser_open("http://gretl.ecn.wfu.edu/apply/");
}

static void login_dialog (login_info *linfo, GtkWidget *parent)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *vbox, *hbox;
    int i;

    login_init(linfo);

    linfo->dlg = gretl_dialog_new(_("gretl: upload"), parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(linfo->dlg));

    hbox = label_hbox(vbox, _("Upload function package"));

    tbl = gtk_table_new(2, 2, FALSE);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 5);

    for (i=0; i<2; i++) {
	char *src = (i == 0)? linfo->login : linfo->pass;
	GtkWidget *entry;

	label = gtk_label_new((i == 0)? _("Login") : _("Password"));
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, i, i+1,
			 GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL,
			 5, 5);

	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 34);
	gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	if (src != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(entry), src);
	}

	if (i == 0) {
	    linfo->login_entry = entry;
	} else {
	    gtk_entry_set_visibility(GTK_ENTRY(entry), FALSE);
	    linfo->pass_entry = entry;
	}
    }

    hbox = label_hbox(vbox, 
		      _("If you don't have a login to the gretl server\n"
			"please see http://gretl.ecn.wfu.edu/apply/.\n"
			"The 'Website' button below should open this page\n"
			"in your web browser."));

    /* control button area */

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(linfo->dlg));

    /* Cancel */
    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), linfo->dlg);

    /* OK */
    button = ok_validate_button(hbox, &linfo->canceled, NULL);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(login_finalize), linfo);
    gtk_widget_grab_default(button);

    /* Website */
    button = gtk_button_new_with_label("Website");
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       button, TRUE);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(web_get_login), NULL);

    gtk_widget_show_all(linfo->dlg);
}

/* Check the package file against gretlfunc.dtd. We call
   this automatically before uploading a package file to
   the server.  The user can also choose to run the test by
   clicking the "Validate" button in the package editor;
   in that case we set the @verbose flag.
*/

static int validate_package_file (const char *fname, int verbose)
{
    const char *gretldir = gretl_home();
    char dtdname[FILENAME_MAX];
    xmlDocPtr doc;
    xmlDtdPtr dtd;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, NULL, &doc, NULL);
    if (err) {
	gui_errmsg(err);
	return 1;
    }

    sprintf(dtdname, "%sfunctions%cgretlfunc.dtd", gretldir, SLASH);
    dtd = xmlParseDTD(NULL, (const xmlChar *) dtdname); 

    if (dtd == NULL) {
	if (verbose) {
	    errbox("Couldn't open DTD to check package");
	} else {
	    fprintf(stderr, "Couldn't open DTD to check package\n");
	}
    } else {
	const char *pkgname = path_last_element(fname);
	xmlValidCtxtPtr cvp = xmlNewValidCtxt();
	PRN *prn = NULL;
	int xerr = 0;

	if (cvp == NULL) {
	    xerr = 1;
	    if (verbose) nomem();
	} else {
	    xerr = bufopen(&prn);
	}

	if (xerr) {
	    xmlFreeDtd(dtd);
	    xmlFreeDoc(doc);
	    return 0;
	}

	cvp->userData = (void *) prn;
	cvp->error    = (xmlValidityErrorFunc) pprintf;
	cvp->warning  = (xmlValidityWarningFunc) pprintf;

	if (!xmlValidateDtd(cvp, doc, dtd)) {
	    const char *buf = gretl_print_get_buffer(prn);

	    errbox(buf);
	    err = 1;
	} else if (verbose) {
	    infobox_printf(_("%s: validated against DTD OK"), pkgname);
	} else {
	    fprintf(stderr, "%s: validated against DTD OK\n", pkgname);
	}

	gretl_print_destroy(prn);
	xmlFreeValidCtxt(cvp);
	xmlFreeDtd(dtd);
    } 

    xmlFreeDoc(doc);

    return err;
}

/* Collect pkg.gfn plus additional package files (PDF doc and/or data
   files) into a temporary dir under the user's dotdir, and make a zip
   archive. If @dest is non-NULL, that's the name of the zipfile to
   build; otherwise if @pzipname is non-NULL the zipfile will be named
   automatically based on the package name, and this name will be
   "returned" in @pzipname.
*/

static int gui_pkg_make_zipfile (function_info *finfo, int pdfdoc,
				 gchar **pzipname, const char *dest)
{
    windata_t *vwin;
    PRN *prn = NULL;
    int err = 0;

    if (pzipname == NULL && dest == NULL) {
	/* we need one or the other */
	return E_DATA;
    }

    /* open printer for recording */
    bufopen(&prn);

    err = package_make_zipfile(finfo->fname,
			       pdfdoc,
			       finfo->datafiles,
			       finfo->n_files,
			       pzipname, dest,
			       OPT_G, prn);
    
    /* show details of operation */
    vwin = view_buffer(prn, 78, 300, _("build zip file"), BUILD_PKG, NULL);
    gtk_window_set_transient_for(GTK_WINDOW(vwin->main),
				 GTK_WINDOW(finfo->dlg));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(vwin->main), TRUE);

    return err;
}

static void do_pkg_upload (function_info *finfo)
{
    const char *fname;
    gchar *buf = NULL;
    char *retbuf = NULL;
    gchar *zipname = NULL;
    login_info linfo;
    GdkDisplay *disp;
    GdkCursor *cursor;
    GdkWindow *w1;
    gint x, y;
    gsize buflen;
    int error_printed = 0;
    int pdfdoc;
    int err;

    err = validate_package_file(finfo->fname, 0);
    if (err) {
	return;
    }

    pdfdoc = pkg_has_pdf_doc(finfo);

    if (pdfdoc || finfo->datafiles != NULL) {
	err = gui_pkg_make_zipfile(finfo, pdfdoc, &zipname, NULL);
	if (err) {
	    /* the error message will have been handled abive */
	    return;
	}
    }

    fname = zipname != NULL ? zipname : finfo->fname;

    login_dialog(&linfo, finfo->dlg);

    if (linfo.canceled) {
	linfo_free(&linfo);
	g_free(zipname);
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

    err = gretl_file_get_contents(fname, &buf, &buflen);

    if (err) {
	error_printed = 1;
    } else {
	err = upload_function_package(linfo.login, linfo.pass, 
				      path_last_element(fname),
				      buf, buflen, &retbuf);
	fprintf(stderr, "upload_function_package: err = %d\n", err);
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

    if (zipname != NULL) {
	/* delete the upload zipfile */
	gretl_remove(zipname);
	g_free(zipname);
    }
    
    g_free(buf);
    free(retbuf);

    linfo_free(&linfo);
}

static int upload_precheck_gfn (const char *fname,
				gchar **zname)
{
    char **datafiles = NULL;
    int pdfdoc = 0;
    int n_files = 0;
    int err;
	    
    err = validate_package_file(fname, 0);
    if (err) {
	return err;
    }
    
    if (package_needs_zipping(fname, &pdfdoc, &datafiles, &n_files)) {
	int resp;

	resp = yes_no_dialog(NULL,
			     _("This package must be uploaded as a zip file.\n"
			       "Try to create the zip file now?"),
			     NULL);
	if (resp == GRETL_YES) {
	    PRN *prn = NULL;

	    if (bufopen(&prn)) {
		err = 1;
	    } else {
		err = package_make_zipfile(fname, pdfdoc,
					   datafiles, n_files,
					   zname, NULL,
					   OPT_G, prn);
		/* show result on error */
		if (err) {
		    view_buffer(prn, 78, 300, _("build zip file"),
				BUILD_PKG, NULL);
		} else {
		    gretl_print_destroy(prn);
		}
	    }
	} else {
	    /* not really an error, but canceled */
	    err = 1;
	}

	strings_array_free(datafiles, n_files);
    }

    return err;
}

/* callback from file selector, so @fname will be a full path */

void upload_specified_package (const char *fname)
{
    const char *realname;
    gchar *zname = NULL;
    gchar *buf = NULL;
    char *retbuf = NULL;
    login_info linfo;
    GdkDisplay *disp;
    GdkCursor *cursor;
    GdkWindow *w1;
    gint x, y;
    gsize buflen;
    int error_printed = 0;
    int err;

    if (has_suffix(fname, ".gfn")) {
	err = upload_precheck_gfn(fname, &zname);
	if (err) {
	    return;
	}
    }

    login_dialog(&linfo, mdata->main);

    if (linfo.canceled) {
	if (zname != NULL) {
	    gretl_remove(zname);
	    g_free(zname);
	}
	linfo_free(&linfo);
	return;
    }

    realname = zname != NULL ? zname : fname;

    /* set waiting cursor */
    disp = gdk_display_get_default();
    w1 = gdk_display_get_window_at_pointer(disp, &x, &y);
    if (w1 != NULL) {
	cursor = gdk_cursor_new(GDK_WATCH);
	gdk_window_set_cursor(w1, cursor);
	gdk_display_sync(disp);
	gdk_cursor_unref(cursor);
    }

    err = gretl_file_get_contents(realname, &buf, &buflen);

    if (err) {
	error_printed = 1;
    } else {
	err = upload_function_package(linfo.login, linfo.pass, 
				      path_last_element(realname),
				      buf, buflen, &retbuf);
	fprintf(stderr, "upload_function_package: err = %d\n", err);
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

    if (zname != NULL) {
	/* we created an on-the-fly temporary zipfile */
	gretl_remove(zname);
	g_free(zname);
    }

    g_free(buf);
    free(retbuf);

    linfo_free(&linfo);
}

/* the basename of a function package file must meet some sanity
   requirements */

static int check_package_filename (const char *fname,
				   int fullpath,
				   GtkWidget *parent)
{
    const char *p = fname;
    int n, err = 0;

    if (fullpath) {
	p = strrchr(fname, SLASH);
	if (p == NULL) {
	    p = fname;
	} else {
	    p++;
	}
    }

    if (fullpath && !has_suffix(p, ".gfn")) {
	/* must have the right suffix */
	err = 1;
    } else {
	n = strlen(p) - (fullpath ? 4 : 0);
	if (n >= FN_NAMELEN) {
	    /* too long */
	    err = 1;
	} else if (gretl_namechar_spn(p) != n) {
	    /* contains funny stuff */
	    err = 1;
	}
    }

    if (err) {
	if (fullpath) {
	    msgbox(_("Invalid package filename: the name must start with a letter,\n"
		     "must be less than 32 characters in length, must include only\n"
		     "ASCII letters, numbers and '_', and must end with \".gfn\"."),
		   GTK_MESSAGE_ERROR, parent);
	} else {
	    msgbox(_("Invalid package name: the name must start with a letter,\n"
		     "must be less than 32 characters in length, and must include\n"
		     "only ASCII letters, numbers and '_'."),
		   GTK_MESSAGE_ERROR, parent);
	}	    
    }

    return err;
}

static int pkg_save_special_functions (function_info *finfo)
{
    const char *key;
    int i, role, err = 0;

    for (i=0; i<N_SPECIALS && !err; i++) {
	role = i + 1;
	key = package_role_get_key(role);
	err = function_set_package_role(finfo->specials[i], 
					finfo->pkg,
					key,
					NULL);
	if (!err && role == UFUN_GUI_MAIN) {
	    /* ensure that the gui-main for a model-window
	       package is set as menu-only */
	    if (finfo->menupath != NULL &&
		strstr(finfo->menupath, "MODELWIN")) {
		finfo->gui_attrs[i] |= UFUN_MENU_ONLY;
	    }
	}
    }

    return err;
}

/* We're saving a previously saved/installed package, and it
   (now) ought to be in its own subdir (PDF doc or data files
   have been specified). We check to see if the gfn file is 
   actually just sitting in /some/path/functions.

   If so, we try to move it into its own subdir and adjust
   everything that depends on its path accordingly.
*/

static int maybe_fix_package_location (function_info *finfo)
{
    const char *pkgname;
    int err = 0;

    pkgname = function_package_get_name(finfo->pkg);

    if (pkg_path_is_toplevel(finfo, pkgname)) {
	char *p, newpath[FILENAME_MAX];

	strcpy(newpath, finfo->fname);
	/* trim off pkgname.gfn */
	p = strrchr(newpath, SLASH);
	*(p+1) = '\0';
	/* append own subdir name */
	strcat(newpath, pkgname);
	/* make/verify the subdir */
	err = gretl_mkdir(newpath);
	if (!err) {
	    /* append pkgname.gfn */
	    strcat(newpath, SLASHSTR);
	    strcat(newpath, pkgname);
	    strcat(newpath, ".gfn");
	    /* and try moving the file */
	    err = gretl_rename(finfo->fname, newpath);
	}
	if (!err) {
	    /* maybe revise "recent" gfn list */
	    delete_from_filelist(FILE_LIST_GFN, finfo->fname);
	    /* update the record in @finfo */
	    g_free(finfo->fname);
	    finfo->fname = g_strdup(newpath);
	    /* and also the in-memory package */
	    function_package_set_properties(finfo->pkg, "fname",
					    newpath, NULL);
	}

	fprintf(stderr, "maybe_fix_package_location: err = %d\n", err);
    } 

    return err;
}

/* Callback from file selector when saving a function package, or
   directly from the package editor if using the package's
   existing filename.
*/

int save_function_package (const char *fname, gpointer p)
{
    function_info *finfo = p;
    int err = 0;

    if (finfo->fname == NULL) {
	/* new save: no filename recorded yet */
	err = check_package_filename(fname, 1, finfo->dlg);
	if (err) {
	    return err;
	}
	finfo->fname = g_strdup(fname);
    }	

    if (finfo->pkg == NULL) {
	/* starting from scratch */
	finfo->pkg = function_package_new(fname, finfo->pubnames, finfo->n_pub,
					  finfo->privnames, finfo->n_priv,
					  &err);
    } else {
	/* revising an existing package */
	err = function_package_connect_funcs(finfo->pkg, finfo->pubnames, finfo->n_pub,
					     finfo->privnames, finfo->n_priv);
	if (err) {
	    fprintf(stderr, "function_package_connect_funcs: err = %d\n", err);
	}
    }

    if (!err) {
	/* we need to do this before setting the "gui-attrs" below */
	pkg_save_special_functions(finfo);
    }

    if (!err) {
	err = function_package_set_properties(finfo->pkg,
					      "author",  finfo->author,
					      "email",   finfo->email,
					      "version", finfo->version,
					      "date",    finfo->date,
					      "description", finfo->pkgdesc,
					      "help", finfo->help,
					      "sample-script", finfo->sample,
					      "data-requirement", finfo->dreq,
					      "min-version", finfo->minver,
					      "menu-attachment", finfo->menupath,
					      "label", finfo->menulabel,
					      "gui-help", finfo->gui_help,
					      "gui-attrs", finfo->gui_attrs,
					      "lives-in-subdir", finfo->uses_subdir,
					      NULL);
	if (err) {
	    fprintf(stderr, "function_package_set_properties: err = %d\n", err);
	}
    }

    if (!err) {
	function_package_set_data_files(finfo->pkg,
					finfo->datafiles,
					finfo->n_files);
    }

    if (!err && finfo->uses_subdir) {
	maybe_fix_package_location(finfo);
    }    

    if (!err) {
	err = function_package_write_file(finfo->pkg);
	if (err) {
	    fprintf(stderr, "function_package_write_file: err = %d\n", err);
	}
    }

    if (err) {
	gui_errmsg(err);
    } else {
	const char *pkgname = function_package_get_name(finfo->pkg);
	gchar *title;

	title = g_strdup_printf("gretl: %s", pkgname);
	gtk_window_set_title(GTK_WINDOW(finfo->dlg), title);
	g_free(title);
	finfo_set_modified(finfo, FALSE);
	gtk_widget_set_sensitive(finfo->validate, TRUE);
	maybe_update_gfn_browser(pkgname,
				 finfo->version,
				 finfo->pkgdesc,
				 finfo->fname,
				 finfo->uses_subdir,
				 pkg_has_pdf_doc(finfo));
	mkfilelist(FILE_LIST_GFN, finfo->fname);
	
	/* revise stored gui package info in accordance with any
	   changes above, as needed */
	gui_function_pkg_revise_status(pkgname,
				       finfo->fname,
				       finfo->menulabel,
				       finfo->menupath,
				       finfo->uses_subdir);
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

    /* write basic package info */
    pprintf(prn, "# author='%s'\n", finfo->author);
    if (finfo->email != NULL && *finfo->email != '\0') {
	pprintf(prn, "# email='%s'\n", finfo->email);
    }
    pprintf(prn, "# version='%s'\n", finfo->version);
    pprintf(prn, "# date='%s'\n", finfo->date);

    /* write private functions, if any */
    for (i=0; i<finfo->n_priv; i++) {
	fun = get_function_from_package(finfo->privnames[i],
					finfo->pkg);
	if (fun != NULL) {
	    pputc(prn, '\n');
	    gretl_function_print_code(fun, tabwidth, prn);
	}
    }

    /* write public functions */
    for (i=0; i<finfo->n_pub; i++) {
	fun = get_function_from_package(finfo->pubnames[i],
					finfo->pkg);
	if (fun != NULL) {
	    pputc(prn, '\n');
	    gretl_function_print_code(fun, tabwidth, prn);
	}
    }

    /* append sample script? */
    if ((finfo->save_flags & APPEND_SAMPLE) &&
	finfo->sample != NULL) {
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

static void maybe_print (PRN *prn, const char *key, 
			 const char *arg)
{
    if (arg != NULL && *arg != '\0') {
	pprintf(prn, "%s = %s\n", key, arg);
    } else {
	pprintf(prn, "%s = \n", key);
    }
}

static int maybe_write_aux_file (function_info *finfo,
				 const char *fname,
				 const char *id,
				 const gchar *content,
				 PRN *prn)
{
    int ret = 0;

    if (content != NULL && *content != '\0') {
	const char *auxname = NULL;
	int flag;

	if (!strcmp(id, "sample-script")) {
	    auxname = finfo->sample_fname;
	    flag = WRITE_SAMPFILE;
	} else if (!strcmp(id, "help")) {
	    auxname = finfo->help_fname;
	    flag = WRITE_HELPFILE;
	} else {
	    auxname = finfo->gui_help_fname;
	    flag = WRITE_GUI_HELP;
	}

	if (auxname != NULL && (finfo->save_flags & flag)) {
	    /* we'll write out the actual file */
	    FILE *fp = NULL;

	    if (strrchr(fname, SLASH)) {
		/* package fname has directory component */
		char *s, tmp[FILENAME_MAX];
		
		strcpy(tmp, fname);
		s = strrchr(tmp, SLASH);
		*(s + 1) = '\0';
		strcat(tmp, auxname);
		fp = gretl_fopen(tmp, "w");
	    } else {
		fp = gretl_fopen(auxname, "w");
	    }
	    
	    if (fp != NULL) {
		fputs(content, fp);
		fputc('\n', fp);
		fclose(fp);
		ret = 1;
	    }
	}

	if (auxname != NULL) {
	    /* record filename in spec file */
	    pputs(prn, auxname);
	}
    }

    return ret;
}

/* Given the in-memory representation of a gfn package, write
   out the corresponding .spec file. Also write out to separate
   files the package's help text and sample script, if available.
*/

int save_function_package_spec (const char *fname, gpointer p)
{
    const char *extra_keys[] = {
	GUI_MAIN,
	"label",
	"menu-attachment",
	BUNDLE_PRINT,
	BUNDLE_PLOT,
	BUNDLE_TEST,
	BUNDLE_FCAST,
	BUNDLE_EXTRA,
	GUI_PRECHECK,
	NULL
    };
    const char *gui_help;
    const char *sample;
    gchar *strval;
    int intval;
    function_info *finfo = p;
    PRN *prn;
    const char *reqstr = NULL;
    int v1, v2, v3;
    int nnp = 0, nmo = 0;
    int i, len;
    int err = 0;

    prn = gretl_print_new_with_filename(fname, &err);
    if (err) {
	file_write_errbox(fname);
	return err;
    }

    maybe_print(prn, "author", finfo->author);
    maybe_print(prn, "email", finfo->email);
    maybe_print(prn, "version", finfo->version);
    maybe_print(prn, "date", finfo->date);
    maybe_print(prn, "description", finfo->pkgdesc);

    get_maj_min_pl(finfo->minver, &v1, &v2, &v3);
    pprintf(prn,"min-version = %d.%d.%d\n", v1, v2, v3); 

    if (finfo->dreq == FN_NEEDS_TS) {
	reqstr = NEEDS_TS;
    } else if (finfo->dreq == FN_NEEDS_QM) {
	reqstr = NEEDS_QM;
    } else if (finfo->dreq == FN_NEEDS_PANEL) {
	reqstr = NEEDS_PANEL;
    } else if (finfo->dreq == FN_NODATA_OK) {
	reqstr = NO_DATA_OK;
    } 

    if (reqstr != NULL) {
	pprintf(prn, "data-requirement = %s\n", reqstr);
    }

    for (i=0; extra_keys[i] != NULL; i++) {
	function_package_get_properties(finfo->pkg, extra_keys[i],
					&strval, NULL);
	if (strval != NULL) {
	    if (*strval != '\0') {
		pprintf(prn, "%s = %s\n", extra_keys[i], strval);
	    }
	    g_free(strval);
	    strval = NULL;
	}
    }

    function_package_get_properties(finfo->pkg, "model-requirement",
				    &intval, NULL);
    if (intval > 0) {
	reqstr = gretl_command_word(intval);
	if (*reqstr != '\0') {
	    pprintf(prn, "model-requirement = %s\n", reqstr);
	}
    }

    /* public interface names */
    pputs(prn, "public = ");
    len = 9;
    for (i=0; i<finfo->n_pub; i++) {
	const char *s = finfo->pubnames[i];
	int n = strlen(s);
	ufunc *fun;

	len += n;
	if (len > 72) {
	    pputs(prn, "\\\n");
	    pprintf(prn, "  %s ", s);
	    len = n + 3;
	} else {
	    pprintf(prn, "%s ", s);
	    len++;
	}
	fun = get_function_from_package(s, finfo->pkg);
	if (user_func_is_noprint(fun)) {
	    nnp++;
	}
	if (user_func_is_menu_only(fun)) {
	    nmo++;
	}	
    }
    pputc(prn, '\n');

    if (nnp > 0) {
	/* no-print interface names */
	pputs(prn, "no-print = ");
	for (i=0; i<finfo->n_pub; i++) {
	    const char *s = finfo->pubnames[i];
	    ufunc *fun = get_function_from_package(s, finfo->pkg);

	    if (user_func_is_noprint(fun)) {
		pprintf(prn, "%s ", s);
	    }
	}
	pputc(prn, '\n');
    }

    if (nmo > 0) {
	/* menu-only interface names */
	pputs(prn, "menu-only = ");
	for (i=0; i<finfo->n_pub; i++) {
	    const char *s = finfo->pubnames[i];
	    ufunc *fun = get_function_from_package(s, finfo->pkg);

	    if (user_func_is_menu_only(fun)) {
		pprintf(prn, "%s ", s);
	    }
	}
	pputc(prn, '\n');
    }

    /* write out help text? */
    pputs(prn, "help = ");
    if (finfo->text != NULL) {
	gchar *help = textview_get_trimmed_text(finfo->text);

	if (help != NULL) {
	    if (!strncmp(help, "pdfdoc", 6)) {
		pputs(prn, help + 7);
	    } else {
		maybe_write_aux_file(finfo, fname, "help",
				     help, prn);
	    }
	    g_free(help);
	}
	pputc(prn, '\n');
    }

    gui_help = (finfo->gui_help != NULL)? finfo->gui_help :
	function_package_get_string(finfo->pkg, "gui-help");

    /* write out GUI-specific help text? */
    if (gui_help != NULL) {
	pputs(prn, "gui-help = ");
	maybe_write_aux_file(finfo, fname, "gui-help",
			     gui_help, prn);
	pputc(prn, '\n');
    }

    sample = (finfo->sample != NULL)? finfo->sample :
	function_package_get_string(finfo->pkg, "sample-script");

    /* write out sample script? */
    pputs(prn, "sample-script = ");
    maybe_write_aux_file(finfo, fname, "sample-script",
			 sample, prn);
    pputc(prn, '\n');

    gretl_print_destroy(prn);

    return 0;
}

int save_function_package_zipfile (const char *fname, gpointer p)
{
    function_info *finfo = p;
    int pdfdoc = pkg_has_pdf_doc(finfo);
    
    gui_pkg_make_zipfile(finfo, pdfdoc, NULL, fname);

    return 0;
}

static int get_lists_from_selector (int **l1, int **l2)
{
    gchar *liststr = get_selector_storelist();
    int err = 0;

    if (liststr == NULL) {
	err = E_DATA;
    } else {
	int *list = gretl_list_from_string(liststr, &err);
    
	if (list != NULL) {
	    if (gretl_list_has_separator(list)) {
		err = gretl_list_split_on_separator(list, l1, l2);
	    } else {
		*l1 = list;
		list = NULL;
	    }
	    free(list);
	}
	g_free(liststr);
    }

    return err;
}

/* called from function selection dialog: a set of functions has been
   selected and now we need to add info on author, version, etc.
*/

void edit_new_function_package (void)
{
    function_info *finfo;
    int *publist = NULL;
    int *privlist = NULL;
    int err = 0;

    err = get_lists_from_selector(&publist, &privlist);
    if (err) {
	return;
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    if (!err) {
	err = finfo_set_function_names(finfo, publist, privlist);
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

/* callback from GUI selector to add/remove functions 
   when editing a package */

void revise_function_package (void *p)
{
    function_info *finfo = p;
    int *pub_ids = NULL;
    int *priv_ids = NULL;
    int changed = 0;
    int err = 0;

    err = get_lists_from_selector(&pub_ids, &priv_ids);
    if (err) {
	return;
    }

#if PKG_DEBUG
    printlist(pub_ids, "new pub_ids");
    printlist(priv_ids, "new priv_ids");
#endif

    fprintf(stderr, "original: finfo->n_pub=%d, finfo->n_priv=%d\n", 
	    finfo->n_pub, finfo->n_priv);

    err = finfo_reset_function_names(finfo, pub_ids, priv_ids,
				     &changed);

    fprintf(stderr, "revised: finfo->n_pub=%d, finfo->n_priv=%d\n", 
	    finfo->n_pub, finfo->n_priv);

    if (!err && changed) {
	depopulate_combo_box(GTK_COMBO_BOX(finfo->codesel));
	func_selector_set_strings(finfo, finfo->codesel);
	verify_selected_specials(finfo);
	if (finfo->pkg != NULL) {
	    /* sync with gretl_func.c */
	    function_package_connect_funcs(finfo->pkg, 
					   finfo->pubnames,
					   finfo->n_pub,
					   finfo->privnames,
					   finfo->n_priv);
	}
	finfo_set_modified(finfo, TRUE);
    }

    free(pub_ids);
    free(priv_ids);
}

static void finfo_set_menuwin (function_info *finfo)
{
    if (finfo->menupath == NULL) {
	finfo->menuwin = NO_WINDOW;
    } else if (!strncmp(finfo->menupath, "MAINWIN", 7)) {
	finfo->menuwin = MAIN_WINDOW;
    } else if (!strncmp(finfo->menupath, "MODELWIN", 8)) {
	finfo->menuwin = MODEL_WINDOW;
    } else {
	finfo->menuwin = NO_WINDOW;
    }
}

static int finfo_set_special_names (function_info *finfo)
{
    const char *key;
    int i, err = 0;

    for (i=0; i<N_SPECIALS && !err; i++) {
	key = package_role_get_key(i+1);
	err = function_package_get_properties(finfo->pkg, key,
					      &finfo->specials[i],
					      NULL);
    }

    return err;
}

static int finfo_set_data_files (function_info *finfo)
{
    char **S;
    int n = 0;

    S = function_package_get_data_files(finfo->pkg, &n);

    if (S != NULL) {
	finfo->datafiles = S;
	finfo->n_files = n;
    }

    return 0;
}

void edit_function_package (const char *fname)
{
    function_info *finfo = NULL;
    int *publist = NULL;
    int *privlist = NULL;
    fnpkg *pkg;
    int err = 0;

    pkg = get_function_package_by_filename(fname, &err);
    if (err) {
	gui_errmsg(err);
	goto bailout;
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    finfo->pkg = pkg;

    err = function_package_get_properties(finfo->pkg,
					  "publist",  &publist,
					  "privlist", &privlist,
					  "author",   &finfo->author,
					  "email",    &finfo->email,
					  "version",  &finfo->version,
					  "date",     &finfo->date,
					  "description", &finfo->pkgdesc,
					  "help", &finfo->help,
					  "sample-script", &finfo->sample,
					  "data-requirement", &finfo->dreq,
					  "min-version", &finfo->minver,
					  "menu-attachment", &finfo->menupath,
					  "label", &finfo->menulabel,
					  "gui-help", &finfo->gui_help,
					  "lives-in-subdir", &finfo->uses_subdir,
					  "gui-attrs", finfo->gui_attrs,
					  NULL);

    if (!err && publist == NULL) {
	err = E_DATA;
    }

    if (!err) {
	err = finfo_set_function_names(finfo, publist, privlist);
    }

    if (!err) {
	err = finfo_set_special_names(finfo);
    }

    if (!err) {
	finfo_set_menuwin(finfo);
    }

    if (!err) {
	err = finfo_set_data_files(finfo);
    }

#if PKG_DEBUG
    printlist(publist, "publist");
    printlist(privlist, "privlist");
#endif

    free(publist);
    free(privlist);

    if (err) {
	fprintf(stderr, "function_package_get_info: failed on %s\n", fname);
	errbox("Couldn't get function package information");
	finfo_free(finfo);
	goto bailout;
    } 

    finfo->fname = g_strdup(fname);

 bailout:

    if (err) {
	delete_from_filelist(FILE_LIST_GFN, fname);
    } else {
	/* record opening */
	gchar *tmp = g_strdup(fname);

	mkfilelist(FILE_LIST_GFN, tmp);
	g_free(tmp);
	
	/* and go for it */
	finfo_dialog(finfo);
    }
}

gboolean edit_specified_package (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");
    gboolean ret = FALSE;
    
    if (fp == NULL) {
	file_read_errbox(fname);
	delete_from_filelist(FILE_LIST_GFN, fname);
    } else {
	fclose(fp);
	edit_function_package(fname);
	ret = TRUE;
    }

    return ret;
}

int no_user_functions_check (GtkWidget *parent)
{
    int err = 0;

    if (n_free_functions() == 0) {
	const gchar *query =
	    N_("No functions are available for packaging at present.\n"
	       "Do you want to write a function now?");
	int resp;

	err = 1;
	resp = yes_no_dialog(_("gretl: function packages"),
			     _(query), parent);
	if (resp == GRETL_YES) {
	    do_new_script(FUNC, NULL);
	}
    } 

    return err;
}

/* called from toolbar.c in response to the "build" option from
   window editing a .spec file */

void build_package_from_spec_file (windata_t *vwin)
{
    char inpname[FILENAME_MAX];
    char gfnname[FILENAME_MAX];
    int resp, err = 0;

    switch_ext(inpname, vwin->fname, "inp");
    err = gretl_test_fopen(inpname, "r");
    if (err) {
	gchar *msg = g_strdup_printf(_("Couldn't open %s"), inpname);
	
	msgbox(msg, GTK_MESSAGE_ERROR, vwin->main);
	g_free(msg);
	return;
    }
    
    switch_ext(gfnname, vwin->fname, "gfn");
    resp = overwrite_gfn_check(gfnname, vwin->main, NULL);

    if (resp == GRETL_YES) {
	int save_batch = gretl_in_batch_mode();
	windata_t *prnwin;
	PRN *prn;

	if (bufopen(&prn)) {
	    return;
	}

	function_package_unload_by_filename(gfnname);

	pprintf(prn, "Found script file '%s'\n", inpname);
	err = execute_script(inpname, NULL, prn, SCRIPT_EXEC | INCLUDE_EXEC,
			     vwin->main);
	if (!err) {
	    err = create_and_write_function_package(gfnname, OPT_G, prn);
	    if (err) {
		pputs(prn, "Failed to produce gfn file\n");
	    } else {
		pprintf(prn, "Wrote '%s'\n", gfnname);
	    }
	}
	gretl_set_batch_mode(save_batch);
	prnwin = view_buffer(prn, 78, 450, _("build gfn file"), BUILD_PKG, NULL);
	gtk_window_set_transient_for(GTK_WINDOW(prnwin->main),
				     GTK_WINDOW(vwin->main));
	gtk_window_set_destroy_with_parent(GTK_WINDOW(prnwin->main), TRUE);
    }
}
