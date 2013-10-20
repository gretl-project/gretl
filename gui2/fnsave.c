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
#include "gretl_xml.h"
#include "libset.h"
#include "dlgutils.h"
#include "datafiles.h"
#include "textbuf.h"
#include "fileselect.h"
#include "gretl_www.h"
#include "winstack.h"
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

#define PKG_DEBUG 0
#define NENTRIES 4

enum {
    NO_WINDOW,
    MAIN_WINDOW,
    MODEL_WINDOW
};

typedef struct function_info_ function_info;
typedef struct login_info_ login_info;

/* package-editing dialog and associated info */

struct function_info_ {
    GtkWidget *dlg;        /* editing dialog box */
    GtkWidget *entries[NENTRIES]; /* entry boxes for author, etc. */
    GtkWidget *text;       /* text box for help */
    GtkWidget *codesel;    /* code-editing selector */
    GtkWidget *save;       /* Save button in dialog */
    GtkWidget *saveas;     /* "Save as" button in dialog */
    GtkWidget *validate;   /* Validate button */
    GtkWidget *popup;      /* popup menu */
    GtkWidget *extra;      /* extra properties child dialog */
    GtkWidget *maintree;   /* main menu selection tree */
    GtkWidget *modeltree;  /* model menu selection tree */
    GtkWidget *currtree;   /* currently displayed menu treeview */
    GtkWidget *alttree;    /* currently undisplayed menu treeview */
    GtkWidget *treewin;    /* scrolled window to hold menu trees */
    windata_t *samplewin;  /* window for editing sample script */
    GList *codewins;       /* list of windows editing function code */
    fnpkg *pkg;            /* pointer to package being edited */
    gchar *fname;          /* package filename */
    gchar *author;         /* package author */
    gchar *version;        /* package version number */
    gchar *date;           /* package last-revised date */
    gchar *pkgdesc;        /* package description */
    gchar *sample;         /* sample script for package */
    gchar *help;           /* package help text */
    char **pubnames;       /* names of public functions */
    char **privnames;      /* names of private functions */
    char **specials;       /* names of special functions */
    int n_pub;             /* number of public functions */
    int n_priv;            /* number of private functions */
    int n_special;         /* (max) number of special functions */
    gboolean uses_subdir;  /* the package has its own subdir (0/1) */
    gchar *menupath;       /* path for menu attachment, if any */
    gchar *menulabel;      /* label for menu attachment, if any */
    int menuwin;           /* code for none/main/model window */
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

static int validate_package_file (const char *fname,
				  int verbose);

static void finfo_set_menuwin (function_info *finfo);

static void edit_code_callback (GtkWidget *w, function_info *finfo);

function_info *finfo_new (void)
{
    function_info *finfo;
    int ns = UFUN_ROLE_MAX - 1;

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) {
	return NULL;
    }

    finfo->specials = strings_array_new(ns);
    if (finfo->specials == NULL) {
	free(finfo);
	return NULL;
    }

    finfo->n_special = ns;

    finfo->pkg = NULL;
    finfo->fname = NULL;
    finfo->author = NULL;
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

    finfo->upload = FALSE;
    finfo->modified = FALSE;

    finfo->active = NULL;
    finfo->samplewin = NULL;
    finfo->codewins = NULL;
    finfo->codesel = NULL;
    finfo->popup = NULL;
    finfo->extra = NULL;

    finfo->help = NULL;
    finfo->pubnames = NULL;
    finfo->privnames = NULL;

    finfo->n_pub = 0;
    finfo->n_priv = 0;

    finfo->dreq = 0;
    finfo->minver = 10804;
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
    g_free(finfo->version);
    g_free(finfo->date);
    g_free(finfo->pkgdesc);
    g_free(finfo->sample);
    g_free(finfo->help);
    g_free(finfo->menupath);
    g_free(finfo->menulabel);

    if (finfo->pubnames != NULL) {
	strings_array_free(finfo->pubnames, finfo->n_pub);
    }
	
    if (finfo->privnames != NULL) {
	strings_array_free(finfo->privnames, finfo->n_priv);
    }

    if (finfo->specials != NULL) {
	strings_array_free(finfo->specials, finfo->n_special);
    }    

    if (finfo->samplewin != NULL) {
	gtk_widget_destroy(finfo->samplewin->main);
    }

    if (finfo->codewins != NULL) {
	g_list_foreach(finfo->codewins, (GFunc) destroy_code_window, NULL);
	g_list_free(finfo->codewins);
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

    pkgname = function_package_get_name(finfo->pkg);

    if (pkgname != NULL && *pkgname != '\0') {
	strcpy(fname, pkgname);
    } else if (finfo->n_pub > 0) {
	strcpy(fname, finfo->pubnames[0]);
    } else {
	strcpy(fname, "pkg");
    }

    if (mode == SAVE_FUNCTIONS_AS) {
	strcat(fname, ".inp");
    } else if (mode == SAVE_GFN_SPEC) {
	strcat(fname, ".spec");
    } else {
	strcat(fname, ".gfn");
    }	
}

/* Check the user-supplied version string for the package: should be
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

/* Callback from the Save button, or "Save as" button, when editing a
   function package.  We first assemble and check the relevant info
   then if the package is new and has not been saved yet (which is
   flagged by finfo->fname being NULL), or if we're called from "Save
   as", we offer a file selector, otherwise we go ahead and save using
   the package's recorded filename.
*/

static int finfo_save (function_info *finfo, int saveas)
{
    const char *missing = "???";
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
	if (*fields[i] == NULL || !strcmp(*fields[i], missing)) {
	    gtk_entry_set_text(GTK_ENTRY(finfo->entries[i]), missing);
	    gtk_editable_select_region(GTK_EDITABLE(finfo->entries[i]), 0, -1);
	    gtk_widget_grab_focus(finfo->entries[i]);
	    return 1;
	} 
    }

    g_free(finfo->help);
    finfo->help = textview_get_trimmed_text(finfo->text);
    if (finfo->help == NULL || !strcmp(finfo->help, missing)) {
	textview_set_text_selected(finfo->text, missing);
	gtk_widget_grab_focus(finfo->text);
	return 1;
    }

    if (finfo->sample == NULL) {
	infobox(_("Please add a sample script for this package"));
	return 1;
    }

    if (check_version_string(finfo->version)) {
	errbox(_("Invalid version string: use numbers and '.' only"));
	return 1;
    }

    if (finfo->fname == NULL || saveas) {
	file_selector_with_parent(SAVE_FUNCTIONS, FSEL_DATA_MISC, 
				  finfo, finfo->dlg);
    } else {
	fprintf(stderr, "Calling save_function_package\n");
	err = save_function_package(finfo->fname, finfo);
    }

    return err;
}

static void finfo_save_callback (GtkWidget *w, function_info *finfo)
{
    int saveas = (w == finfo->saveas);

    finfo_save(finfo, saveas);
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

static int extract_funcname (char *name, const char *s)
{
    char word[FN_NAMELEN];
    int n, type, err = 0;

    s += strspn(s, " ");
    n = strcspn(s, " (");

    if (n == 0 || n > FN_NAMELEN - 1) {
	return E_DATA;
    }

    *word = '\0';
    strncat(word, s, n);

    if (!strcmp(word, "void")) {
	type = GRETL_TYPE_VOID;
    } else {
	type = gretl_type_from_string(word);
    }

    if (type == 0) {
	/* old-style */
	strcpy(name, word);
    } else if (!ok_function_return_type(type)) {
	err = E_DATA;
    } else {
	s += n;
	s += strspn(s, " ");
	n = strcspn(s, " (");
	if (n == 0 || n > FN_NAMELEN - 1) {
	    err = E_DATA;
	} else {
	    *name = '\0';
	    strncat(name, s, n);
	}
    }

    return err;
}

static int pretest_funcname (char *buf, const char *origname)
{
    char *s, line[MAXLINE], newname[FN_NAMELEN];
    int err = 0;

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && !err) {
	s = line + strspn(line, " \t");
	if (!strncmp(s, "function ", 9)) {
	    if (extract_funcname(newname, s + 9)) {
		err = E_DATA;
	    } else if (strcmp(newname, origname)) {
		gretl_errmsg_set(_("You can't change the name of a function here"));
		err = E_DATA;
	    }
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
	err = execute_script(NULL, text, NULL, INCLUDE_EXEC);
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
	    errbox(_("Can't find the function '%s'"), funname);
	}
    }

    if (bufopen(&prn)) {
	return;
    }    

    gretl_function_print_code(fun, prn);

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
	finfo->sample = gretl_strdup(text);
	g_free(text);
	mark_vwin_content_saved(vwin);
	finfo_set_modified(finfo, TRUE);
    }
}

static void nullify_sample_window (GtkWidget *w, function_info *finfo)
{
    finfo->samplewin = NULL;
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

static void add_remove_callback (GtkWidget *w, function_info *finfo)
{
    add_remove_functions_dialog(finfo->pubnames, finfo->n_pub,
				finfo->privnames, finfo->n_priv,
				finfo->pkg, finfo);
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

static void gfn_to_spec_callback (GtkWidget *w, function_info *finfo)
{
    file_selector_with_parent(SAVE_GFN_SPEC, FSEL_DATA_MISC, 
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
				       const int *publist,
				       const int *privlist,
				       int *changed)
{
    char **pubnames = NULL;
    char **privnames = NULL;
    int npub = publist[0];
    int npriv = (privlist == NULL)? 0 : privlist[0];
    int err = 0;

    *changed = 0;

    pubnames = get_function_names(publist, &err);

    if (!err && npriv > 0) {
	privnames = get_function_names(privlist, &err);
    }

    if (!err) {
	if (npub != finfo->n_pub || npriv != finfo->n_priv) {
	    /* we know that something has changed */
	    *changed = 1;
	} else {
	    /* we'll have to check for any changes */
	    *changed = strings_array_cmp(pubnames, finfo->pubnames, npub);
	    if (!*changed && npriv > 0) {
		*changed = strings_array_cmp(privnames, finfo->privnames,
					     npriv);
	    }
	}
    }

    if (!*changed) {
	strings_array_free(pubnames, npub);
	strings_array_free(privnames, npriv);
    } else if (!err) {
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

    get_maj_min_pl(finfo->minver, &maj, &min, &pl);

    tmp = gtk_label_new(_("Minimum gretl version"));
    gtk_misc_set_alignment(GTK_MISC(tmp), 1.0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
    gtk_widget_show(tmp);

    hbox = gtk_hbox_new(FALSE, 0);

    spin = gtk_spin_button_new_with_range(1, 9, 1);
    if (maj > 1) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), (double) maj);
    }
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 2);
    g_object_set_data(G_OBJECT(spin), "level", GINT_TO_POINTER(1));
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(adjust_minver), finfo);
    tmp = gtk_label_new(".");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 2);

    spin = gtk_spin_button_new_with_range(0, 99, 1);
    if (min > 0) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), (double) min);
    }
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 2);
    g_object_set_data(G_OBJECT(spin), "level", GINT_TO_POINTER(2));
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(adjust_minver), finfo);
    tmp = gtk_label_new(".");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 2);

    spin = gtk_spin_button_new_with_range(0, 99, 1);
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
	int resp = yes_no_dialog("gretl", _("Save changes?"), 1);

	if (resp == GRETL_CANCEL) {
	    return TRUE;
	} else if (resp == GRETL_YES) {
	    return finfo_save(finfo, 0);
	}
    }

    return FALSE;
}

static void check_pkg_callback (GtkWidget *widget, function_info *finfo)
{
    validate_package_file(finfo->fname, 1);
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
    gtk_box_pack_start(GTK_BOX(holder), hbox, FALSE, FALSE, 5);
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), finfo->menuwin);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(switch_menu_view), finfo);
}

static int process_menu_attachment (function_info *finfo,
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
	    g_free(finfo->menulabel);
	    finfo->menulabel = NULL;
	    changed = 1;
	}
    } else if (finfo->menulabel == NULL || strcmp(finfo->menulabel, label)) {
	g_free(finfo->menulabel);
	finfo->menulabel = label;
	label = NULL;
	changed = 1;
    }

    g_free(label);

    if (!gtk_widget_is_sensitive(view)) {
	/* no menu attachment at present */
	if (finfo->menupath != NULL) {
	    g_free(finfo->menupath);
	    finfo->menupath = NULL;
	    changed = 1;
	}
	finfo->menuwin = NO_WINDOW;
	goto finish;
    }	

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));

    if (!gtk_tree_selection_get_selected(selection, &model, &iter)) {
	if (finfo->menupath != NULL) {
	    g_free(finfo->menupath);
	    finfo->menupath = NULL;
	    changed = 1;
	}
    } else {
	gchar *newpath = NULL;

	gtk_tree_model_get(model, &iter, 1, &newpath, -1);
	if (finfo->menupath == NULL || strcmp(finfo->menupath, newpath)) {
	    g_free(finfo->menupath);
	    finfo->menupath = newpath;
	    newpath = NULL;
	    changed = 1;
	}
	g_free(newpath);
    }

    finfo_set_menuwin(finfo);

 finish:

    if (finfo->menupath != NULL && *finfo->menupath != '\0' &&
	(finfo->menulabel == NULL || *finfo->menulabel == '\0')) {
	warnbox(_("To create a menu attachment, you must supply a label."));
	*focus_label = 1;
    } else if (changed) {
	infobox(_("To update the menu attachment, you should\n"
		  "(a) save this package, and (b) restart gretl."));
    }

    return changed;
}

static int process_special_functions (function_info *finfo)
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

    for (i=0; i<finfo->n_special && !err; i++) {
	if (gtk_widget_is_sensitive(c_array[i])) {
	    int newnull = 0, oldnull = 0, changed = 0;

	    /* retrieve and check the selected name */
	    newfun = combo_box_get_active_text(GTK_COMBO_BOX(c_array[i]));
	    newnull = (newfun == NULL || *newfun == '\0' || 
		       !strcmp(newfun, "none"));

	    /* retrieve and check what was there before */
	    oldfun = finfo->specials[i];
	    oldnull = (oldfun == NULL || *oldfun == '\0' || 
		       !strcmp(oldfun, "none"));

	    if (oldnull && !newnull) {
		changed = 1;
	    } else if (!oldnull && newnull) {
		changed = 1;
	    } else if (!oldnull && !newnull) {
		changed = strcmp(newfun, oldfun);
	    } 
	    if (changed) {
		free(finfo->specials[i]);
		finfo->specials[i] = gretl_strdup(newfun);
		n_changed++;
	    }
	    g_free(newfun);
	}
    }

    return n_changed;
}

static void extra_properties_callback (GtkWidget *w, function_info *finfo)
{
    GtkWidget *notebook;
    gint page;
    int changed;

    notebook = g_object_get_data(G_OBJECT(finfo->extra), "book");
    page = gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));
    
    if (page == 0) {
	changed = process_special_functions(finfo);
    } else {
	int focus_label = 0;

	changed = process_menu_attachment(finfo, &focus_label);
	if (focus_label) {
	    GtkWidget *w = g_object_get_data(G_OBJECT(finfo->extra), "label-entry");

	    gtk_widget_grab_focus(w);
	}
    }

    if (changed) {
	finfo_set_modified(finfo, TRUE);
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

    s0 = combo_box_get_active_text(this);
    c_array = g_object_get_data(G_OBJECT(finfo->extra), "combo-array");

    for (i=0; i<finfo->n_special && !dup; i++) {
	other = c_array[i];
	if (other != GTK_WIDGET(this)) {
	    si = combo_box_get_active_text(GTK_COMBO_BOX(other));
	    if (!strcmp(si, s0)) {
		/* switch to "none" */
		gtk_combo_box_set_active(GTK_COMBO_BOX(other), 0);
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
    } else {
	show_gui_help(MENU_ATTACH);
    }
}

static void unref_trees (GtkWidget *w, function_info *finfo)
{
    g_object_unref(G_OBJECT(finfo->maintree));
    g_object_unref(G_OBJECT(finfo->modeltree));
}

#define must_be_private(r) (r == UFUN_GUI_PRECHECK)
#define must_be_public(r) (r != UFUN_GUI_PRECHECK)

/* The following function supports two "notebook" tabs: one allows the
   user to select functions in the package for the various "special"
   package roles (e.g. gui-main, bundle-print). The other allows for
   selection of a menu attachment point and GUI label.
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
    int nfuns, i, j;

    if (finfo->extra != NULL) {
	gtk_window_present(GTK_WINDOW(finfo->extra));
	return;
    }

    finfo->maintree = finfo->modeltree = NULL;
    finfo->currtree = finfo->alttree = NULL;
    finfo->treewin = NULL;

    dlg = gretl_dialog_new(_("gretl: extra properties"), finfo->dlg,
			   GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);
    finfo->extra = dlg;
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

    table = gtk_table_new(finfo->n_special, 2, TRUE);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);

    nfuns = finfo->n_priv + finfo->n_pub;
    combo_array = g_malloc(finfo->n_special * sizeof *combo_array);
    g_object_set_data_full(G_OBJECT(dlg), "combo-array", 
			   combo_array, g_free);

    /* For each "special" function role, test the functions
       in finfo->pkg to see if they qualify as candidates for
       that role; if so, add them to the combo selector.
    */

    for (i=0; i<finfo->n_special; i++) {
	int err, n_cands = 0;
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
	    /* test for validity in this role */
	    err = function_set_package_role(funname, finfo->pkg,
					    key, OPT_T, NULL);
	    if (!err) {
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
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), table, FALSE, FALSE, 5);

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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Apply button */
    tmp = apply_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(extra_properties_callback), finfo);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* Close button */
    tmp = close_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
                     G_CALLBACK(delete_widget), dlg);

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
    int i;

    finfo->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_default_size(GTK_WINDOW(finfo->dlg), 640, 500);

    if (finfo->fname != NULL) {
	gchar *title = title_from_filename(finfo->fname, TRUE);

	gtk_window_set_title(GTK_WINDOW(finfo->dlg), title);
	g_free(title);
    } else {
	gtk_window_set_title(GTK_WINDOW(finfo->dlg), 
			     _("gretl: function package editor"));
    } 

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

    hbox = label_hbox(vbox, _("Help text:"));

    finfo->text = editable_text_box(&hbuf);
    text_table_setup(vbox, finfo->text);

    if (finfo->help != NULL) {
	textview_set_text(finfo->text, finfo->help);
    }
    g_signal_connect(G_OBJECT(hbuf), "changed", 
		     G_CALLBACK(pkg_changed), finfo);

    /* edit code button, possibly with function selector */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    button = gtk_button_new_with_label(_("Edit function code"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_code_callback), finfo);

    finfo->codesel = active_func_selector(finfo);
    gtk_box_pack_start(GTK_BOX(hbox), finfo->codesel, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(finfo->codesel), "changed",
		     G_CALLBACK(update_active_func), finfo);

    update_active_func(NULL, finfo);

    /* save package as script button */
    button = gtk_button_new_with_label(_("Save as script"));
    gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(gfn_to_script_callback), finfo);

    /* edit sample script button */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    button = gtk_button_new_with_label(_("Edit sample script"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(edit_sample_callback), finfo);

    /* add/remove functions button */
    button = gtk_button_new_with_label(_("Add/remove functions"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(add_remove_callback), finfo);

    /* extra package properties button */
    button = gtk_button_new_with_label(_("Extra properties"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(extra_properties_dialog), finfo);

    /* write spec file button */
    button = gtk_button_new_with_label(_("Write spec file"));
    gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(gfn_to_spec_callback), finfo);

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

    /* "Validate" button */
    button = gtk_button_new_with_label(_("Validate"));
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(check_pkg_callback), finfo);
    gtk_widget_set_sensitive(button, finfo->fname != NULL);
    finfo->validate = button;

    /* "Save as" button */
    button = gtk_button_new_from_stock(GTK_STOCK_SAVE_AS);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(finfo_save_callback), finfo);
    gtk_widget_set_sensitive(button, finfo->fname != NULL);
    finfo->saveas = button;

    /* Save button */
    button = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(finfo_save_callback), finfo);
    finfo->save = button;

    /* Close button */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(delete_pkg_editor), finfo);

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
	    infobox(_("%s: validated against DTD OK"), pkgname);
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

static void do_pkg_upload (function_info *finfo, const char *fname)
{
    gchar *buf = NULL;
    char *retbuf = NULL;
    login_info linfo;
    GdkDisplay *disp;
    GdkCursor *cursor;
    GdkWindow *w1;
    gint x, y;
    int error_printed = 0;
    int err;

    err = validate_package_file(fname, 0);
    if (err) {
	return;
    }

    login_dialog(&linfo, finfo->dlg);

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

    err = gretl_file_get_contents(fname, &buf, NULL);

    if (err) {
	error_printed = 1;
    } else {
	err = upload_function_package(linfo.login, linfo.pass, 
				      path_last_element(fname),
				      buf, &retbuf);
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

    g_free(buf);
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
	errbox(_("Invalid package filename: the name must start with a letter,\n"
		 "must be less than 32 characters in length, must include only\n"
		 "ASCII letters, numbers and '_', and must end with \".gfn\""));
    }

    return err;
}

static int dont_overwrite_pkg (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");
    int ret = 0;

    if (fp != NULL) {
	int resp = yes_no_dialog("gretl", _("OK to overwrite?"), 0);

	ret = (resp == GRETL_NO);
	fclose(fp);
    }

    return ret;
}

static int pkg_save_special_functions (function_info *finfo)
{
    const char *key;
    int i, err = 0;

    for (i=0; i<finfo->n_special && !err; i++) {
	key = package_role_get_key(i+1);
	err = function_set_package_role(finfo->specials[i], 
					finfo->pkg,
					key,
					OPT_NONE,
					NULL);
    }

    return err;
}

/* callback from file selector when saving a function package, or
   directly from the package editor if using the package's
   existing filename */

int save_function_package (const char *fname, gpointer p)
{
    function_info *finfo = p;
    int err = 0;

    if (finfo->fname == NULL || strcmp(finfo->fname, fname)) {
	/* new or save as */
	if (dont_overwrite_pkg(fname)) {
	    return 1;
	}
	err = check_package_filename(finfo, fname);
	if (err) {
	    return err;
	}
    }	

    if (finfo->fname == NULL) {
	/* no filename recorded yet */
	finfo->fname = g_strdup(fname);
    } else if (strcmp(finfo->fname, fname)) {
	/* doing "Save as" */
	function_package_unload_by_filename(finfo->fname);
	free(finfo->fname);
	finfo->fname = g_strdup(fname);
	finfo->pkg = NULL;
    }

    if (finfo->pkg == NULL) {
	/* starting from scratch, or save as */
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
	err = function_package_set_properties(finfo->pkg,
					      "author",  finfo->author,
					      "version", finfo->version,
					      "date",    finfo->date,
					      "description", finfo->pkgdesc,
					      "help", finfo->help,
					      "sample-script", finfo->sample,
					      "data-requirement", finfo->dreq,
					      "min-version", finfo->minver,
					      "menu-attachment", finfo->menupath,
					      "label", finfo->menulabel,
					      NULL);
	if (err) {
	    fprintf(stderr, "function_package_set_properties: err = %d\n", err);
	}
    }

    if (!err) {
	pkg_save_special_functions(finfo);
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
	gtk_widget_set_sensitive(finfo->saveas, TRUE);
	gtk_widget_set_sensitive(finfo->validate, TRUE);
	maybe_update_func_files_window(EDIT_FN_PKG);
	if (finfo->upload) {
	    do_pkg_upload(finfo, fname);
	}
	/* revise packages.xml in accordance with any
	   changes above */
	revise_package_status(pkgname, 
			      finfo->menulabel,
			      finfo->menupath,
			      finfo->uses_subdir,
			      TRUE, FALSE);
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
	fun = get_function_from_package(finfo->privnames[i],
					finfo->pkg);
	if (fun != NULL) {
	    pputc(prn, '\n');
	    gretl_function_print_code(fun, prn);
	}
    }

    for (i=0; i<finfo->n_pub; i++) {
	fun = get_function_from_package(finfo->pubnames[i],
					finfo->pkg);
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

static void maybe_print (PRN *prn, const char *key, 
			       const char *arg)
{
    if (arg != NULL) {
	pprintf(prn, "%s = %s\n", key, arg);
    } else {
	pprintf(prn, "%s = \n", key);
    }
}

static const char *aux_basename (const char *fname)
{
    const char *p = strrchr(fname, SLASH);

    if (p != NULL) {
	/* take last part of filename */
	return p + 1;
    } else {
	return fname;
    }
}

static int maybe_write_aux_file (const char *fname, 
				 const gchar *content,
				 const char *sfx,
				 PRN *prn)
{
    int ret = 0;

    if (content != NULL && *content != '\0') {
	char *s, tmp[FILENAME_MAX];
	FILE *fp;

	strcpy(tmp, fname);
	s = strstr(tmp, ".spec");
	if (s != NULL) {
	    *s = '\0';
	    strncat(s, sfx, strlen(sfx));
	    fp = gretl_fopen(tmp, "w");
	    if (fp != NULL) {
		fputs(content, fp);
		fclose(fp);
		pputs(prn, aux_basename(tmp));
		ret = 1;
	    }
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
    gchar *strval;
    int intval;
    function_info *finfo = p;
    PRN *prn;
    const char *reqstr = NULL;
    int v1, v2, v3;
    int i, len, nnp = 0;
    int err = 0;

    prn = gretl_print_new_with_filename(fname, &err);
    if (err) {
	file_write_errbox(fname);
	return err;
    }

    maybe_print(prn, "author", finfo->author);
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

    /* write out help text file and/or sample script? */

    pputs(prn, "help = ");
    if (finfo->text != NULL) {
	gchar *help = textview_get_trimmed_text(finfo->text);

	if (help != NULL) {
	    if (!strncmp(help, "pdfdoc:", 7)) {
		pputs(prn, help + 7);
	    } else {
		maybe_write_aux_file(fname, help, "_help.txt", prn);
	    }
	    g_free(help);
	}
	pputc(prn, '\n');
    }

    pputs(prn, "sample-script = ");
    maybe_write_aux_file(fname, finfo->sample, "_sample.inp", prn);
    pputc(prn, '\n');

    gretl_print_destroy(prn);

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
    int *publist = NULL;
    int *privlist = NULL;
    int changed = 0;
    int err = 0;

    err = get_lists_from_selector(&publist, &privlist);
    if (err) {
	return;
    }

#if PKG_DEBUG
    printlist(publist, "new publist");
    printlist(privlist, "new privlist");
#endif

    err = finfo_reset_function_names(finfo, publist, privlist,
				     &changed);

    if (!err && changed) {
	depopulate_combo_box(GTK_COMBO_BOX(finfo->codesel));
	func_selector_set_strings(finfo, finfo->codesel);
	finfo_set_modified(finfo, TRUE);
    }

    fprintf(stderr, "finfo->n_pub=%d, finfo->n_priv=%d\n", 
	    finfo->n_pub, finfo->n_priv);

    free(publist);
    free(privlist);
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

    for (i=0; i<finfo->n_special && !err; i++) {
	key = package_role_get_key(i+1);
	err = function_package_get_properties(finfo->pkg, key,
					      &finfo->specials[i],
					      NULL);
    }

    return err;
}

void edit_function_package (const char *fname)
{
    function_info *finfo;
    int *publist = NULL;
    int *privlist = NULL;
    fnpkg *pkg;
    const char *p;
    int err = 0;

    pkg = get_function_package_by_filename(fname, &err);
    if (err) {
	gui_errmsg(err);
	return;
    }

    finfo = finfo_new();
    if (finfo == NULL) {
	return;
    }

    finfo->pkg = pkg;

    err = function_package_get_properties(finfo->pkg,
					  "publist",  &publist,
					  "privlist", &privlist,
					  "author",   &finfo->author,
					  "version",  &finfo->version,
					  "date",     &finfo->date,
					  "description", &finfo->pkgdesc,
					  "help", &finfo->help,
					  "sample-script", &finfo->sample,
					  "data-requirement", &finfo->dreq,
					  "min-version", &finfo->minver,
					  "menu-attachment", &finfo->menupath,
					  "label", &finfo->menulabel,
					  "lives-in-subdir", &finfo->uses_subdir,
					  NULL);

    if (!err && publist != NULL) {
	err = finfo_set_function_names(finfo, publist, privlist);
    }

    if (!err) {
	err = finfo_set_special_names(finfo);
    }

    if (!err) {
	finfo_set_menuwin(finfo);
    }

#if PKG_DEBUG
    printlist(publist, "publist");
    printlist(privlist, "privlist");
#endif

    free(publist);
    free(privlist);

    if (err || publist == NULL) {
	fprintf(stderr, "function_package_get_info: failed on %s\n", fname);
	errbox("Couldn't get function package information");
	finfo_free(finfo);
	return;
    } 

    /* if the user has a package-list window open, we may need to
       sync (check the "loaded" box for this package if it wasn't
       already set) 
    */
    maybe_update_func_files_window(CALL_FN_PKG);

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
	edit_function_package(fname);
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
