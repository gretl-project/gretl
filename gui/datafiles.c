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

/* datafiles.c : for gretl */

#define COLL_DEBUG 0
#define GFN_DEBUG 0

#include "gretl.h"
#include "gui_utils.h"
#include "datafiles.h"
#include "database.h"
#include "filelists.h"
#include "gretl_www.h"
#include "menustate.h"
#include "fnsave.h"
#include "fncall.h"
#include "treeutils.h"
#include "selector.h"
#include "toolbar.h"
#include "winstack.h"
#include "fileselect.h"
#include "dlgutils.h"
#include "textbuf.h"

#include "gretl_xml.h"
#include "gretl_func.h"
#include "addons_utils.h"

#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

static GtkWidget *files_vbox (windata_t *vwin);
static GtkWidget *files_notebook (windata_t *vwin, int role);
static int populate_notebook_filelists (windata_t *vwin,
					GtkWidget *notebook,
					int role);
static gint populate_gfn_list (windata_t *vwin);

typedef struct _file_collection file_collection;

struct _file_collection {
    char *path;
    char *descfile;
    char *title;
    int type;
    GtkWidget *listbox;
};

enum {
    COLL_DATA,
    COLL_PS,
    COLL_MAX
};

enum {
    PKG_ATTR_RES = 1 << 0,
    PKG_ATTR_DOC = 1 << 1
};

#define GFN_DIRNAME_COL 5
#define GFN_FLAGS_COL 6

#define REMOTE_ACTION(c) (c == REMOTE_DB || \
                          c == REMOTE_FUNC_FILES || \
                          c == REMOTE_DATA_PKGS || \
	                  c == DBNOMICS_DB || \
	                  c == DBNOMICS_SERIES)

#define DBNOMICS_ACTION(c) (c == DBNOMICS_DB || c == DBNOMICS_SERIES)

static GList *collections[COLL_MAX];
static gboolean collections_built;

static int role_to_index (int role)
{
    if (role == TEXTBOOK_DATA) {
	return COLL_DATA;
    } else if (role == PS_FILES) {
	return COLL_PS;
    } else {
	return -1;
    }
}

static GList *collections_for_role (int role)
{
    int i = role_to_index(role);

    if (i >= 0 && collections[i] != NULL) {
	return g_list_first(collections[i]);
    } else {
	return NULL;
    }
}

static void
read_fn_files_in_dir (GDir *dir, const char *path,
		      GtkListStore *store,
		      GtkTreeIter *iter,
		      int *nfn);

static char *full_path (char *s1, const char *s2)
{
    static char fpath[FILENAME_MAX];
    int n = strlen(s1);

    if (s1[n-1] == '.') {
	s1[n-1] = '\0';
	n--;
    }

    if (IS_SLASH(s1[n-1])) {
	sprintf(fpath, "%s%s", s1, s2);
    } else {
	sprintf(fpath, "%s%c%s", s1, SLASH, s2);
    }

#if COLL_DEBUG > 1
    fprintf(stderr, "full_path: got '%s' from '%s' + '%s'\n",
	    fpath, s1, s2);
#endif

    return fpath;
}

/* check for a few known, older, file collections whose
   descriptions files do not conform to the now-standard
   pattern
*/

static int is_oldstyle_collection (file_collection *coll, int *err)
{
    const file_collection std_data[] = {
	{ "wooldridge", "jw_descriptions", "Wooldridge", COLL_DATA, NULL },
	{ "gujarati", "dg_descriptions", "Gujarati", COLL_DATA, NULL },
	{ "pwt56", "descriptions", "PWT 56", COLL_DATA, NULL }
    };
    const file_collection std_ps = {
	"pwt56", "ps_descriptions", "PWT 56", COLL_PS, NULL
    };
    int i;

    for (i=0; i<3; i++) {
	if (strstr(coll->path, std_data[i].path) &&
	    !strcmp(coll->descfile, std_data[i].descfile)) {
	    coll->title = gretl_strdup(std_data[i].title);
	    if (coll->title == NULL) {
		*err = E_ALLOC;
	    } else {
		coll->type = COLL_DATA;
	    }
	    return 1;
	}
    }

    if (strstr(coll->path, std_ps.path) &&
	!strcmp(coll->descfile, std_ps.descfile)) {
	coll->title = gretl_strdup(std_ps.title);
	if (coll->title == NULL) {
	    *err = E_ALLOC;
	} else {
	    coll->type = COLL_PS;
	}
	return 1;
    }

    return 0;
}

/* return non-zero only on fatal error */

static int get_title_from_descfile (file_collection *coll)
{
    char line[64], title[24];
    char *test;
    FILE *fp;
    int err = 0;

    test = full_path(coll->path, coll->descfile);
    fp = gretl_fopen(test, "rb");

    if (fp != NULL && fgets(line, sizeof line, fp) != NULL) {
	gretl_strstrip(line);
	if (sscanf(line, "# %23[^:]", title) == 1) {
	    coll->title = gretl_strdup(title);
	    if (coll->title == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (fp != NULL) {
	fclose(fp);
    }

    return err;
}

static void free_file_collection (file_collection *coll)
{
    free(coll->path);
    free(coll->descfile);
    free(coll->title);
    free(coll);
}

static file_collection *file_collection_new (const char *path,
					     const char *descfile,
					     int *err)
{
    file_collection *coll = malloc(sizeof *coll);

    if (coll == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    coll->type = -1;
    coll->title = NULL;
    coll->path = gretl_strdup(path);
    coll->descfile = gretl_strdup(descfile);

    if (coll->path == NULL || coll->descfile == NULL) {
	*err = E_ALLOC;
    } else {
	int os = is_oldstyle_collection(coll, err);

	if (!*err && !os) {
	    if (strstr(coll->descfile, "ps_")) {
		coll->type = COLL_PS;
	    } else {
		coll->type = COLL_DATA;
	    }
	    *err = get_title_from_descfile(coll);
	}
    }

    if (*err || coll->title == NULL) {
	free_file_collection(coll);
	coll = NULL;
    }

    return coll;
}

static int compare_colls (const void *a, const void *b)
{
    const file_collection *ca = a;
    const file_collection *cb = b;

    if (!strcmp(ca->title, "Gretl")) {
	return -1;
    } else if (!strcmp(cb->title, "Gretl")) {
	return 1;
    } else {
	return strcmp(ca->title, cb->title);
    }
}

static void push_collection (file_collection *collection)
{
    int i = collection->type;

    collections[i] = g_list_append(collections[i], collection);
}

void destroy_file_collections (void)
{
    GList *L;
    int i;

    for (i=0; i<COLL_MAX; i++) {
	if (collections[i] != NULL) {
	    L = g_list_first(collections[i]);
	    while (L) {
		free_file_collection(L->data);
		L = L->next;
	    }
	    g_list_free(collections[i]);
	    collections[i] = NULL;
	}
    }

    collections_built = FALSE;
}

static void sort_files_stack (int role)
{
    int i = role_to_index(role);

    if (i >= 0 && collections[i] != NULL) {
	collections[i] = g_list_sort(collections[i], compare_colls);
    }
}

/* Returns the number of file collections found and pushed;
   writes non-zero to @err if something show-stopping
   occurs
*/

static int get_file_collections_from_dir (const char *path, GDir *dir,
					  int *err)
{
    file_collection *coll;
    const gchar *dname;
    int n = 0;

    while (!*err && (dname = g_dir_read_name(dir))) {
	/* we're looking for a filename that ends with "descriptions" */
	if (strstr(dname, "descriptions")) {
	    size_t len = strlen(dname);

#if COLL_DEBUG
	    fprintf(stderr, "   %s: looking at '%s'\n", path, dname);
#endif
	    if (!strcmp(dname + len - 12, "descriptions")) {
		coll = file_collection_new(path, dname, err);
		if (coll != NULL) {
		    push_collection(coll);
		    n++;
		}
	    }
	}
    }

    return n;
}

static int dont_go_there (const char *s)
{
    int ret = 0;

    if (!strcmp(s, "..") || strstr(s, ".inp") || strstr(s, ".gdt") ||
	strstr(s, ".gretl") || strstr(s, ".hdr")) {
	ret = 1;
    }

    return ret;
}

/* Returns the number of collections found; @err is set to
   non-zero only if something show-stopping occurs
*/

static int seek_file_collections (const char *basedir,
				  SearchType stype,
				  int *err)
{
    gchar *path = NULL;
    GDir *topdir;
    const char *dname;
    int n_coll = 0;

#if COLL_DEBUG
    fprintf(stderr, "*** seek_file_collections: basedir='%s', type=%d\n",
	    basedir, stype);
#endif

    if (stype == DATA_SEARCH) {
	path = g_strdup_printf("%sdata", basedir);
    } else if (stype == SCRIPT_SEARCH) {
	path = g_strdup_printf("%sscripts", basedir);
    } else {
	/* USER_SEARCH */
	path = g_strdup(basedir);
	trim_slash(path);
    }

    topdir = gretl_opendir(path);
    if (topdir == NULL) {
	g_free(path);
	return 0;
    }

#if COLL_DEBUG
    fprintf(stderr, "*** seek_file_collections: path='%s'\n", path);
#endif

    while (!*err && (dname = g_dir_read_name(topdir))) {
	if (!dont_go_there(dname)) {
	    gchar *subpath;
	    GDir *subdir;

#if COLL_DEBUG > 1
	    fprintf(stderr, " dname = '%s'\n", dname);
#endif
	    subpath = g_build_path("/", path, dname, NULL);
	    subdir = gretl_opendir(subpath);
	    if (subdir != NULL) {
#if COLL_DEBUG
		fprintf(stderr, " trying in subdir '%s'\n", subpath);
#endif
		n_coll += get_file_collections_from_dir(subpath, subdir, err);
#if COLL_DEBUG
		if (*err) {
		    fprintf(stderr, " result: err = %d\n", *err);
		}
#endif
		g_dir_close(subdir);
	    }
	    g_free(subpath);
	}
    }

    g_dir_close(topdir);
    g_free(path);

#if COLL_DEBUG
    fprintf(stderr, "*** found %d collections\n", n_coll);
#endif

    return n_coll;
}

#if COLL_DEBUG

static void print_collection (const file_collection *coll)
{
    fprintf(stderr, "path = '%s'\n", coll->path);
    fprintf(stderr, "descfile = '%s'\n", coll->descfile);
    if (coll->title != NULL && *coll->title != '\0') {
	fprintf(stderr, "title = '%s'\n", coll->title);
    }
}

static void print_collections (int role)
{
    GList *L;
    int i;

    if (role == TEXTBOOK_DATA) {
	fputs("\n*** Data collections:\n", stderr);
    } else {
	fputs("\n*** Script collections:\n", stderr);
    }

    for (i=0; i<COLL_MAX; i++) {
	L = g_list_first(collections[i]);
	while (L) {
	    print_collection(L->data);
	    L = L->next;
	}
    }
}

#endif /* COLL_DEBUG */

static int build_file_collections (int role)
{
    static int err;

    if (!collections_built && !err) {
	const char *wd;
	int derr[3] = {0};
	int serr[3] = {0};
	int uerr[2] = {0};
	int nd = 0;
	int ns = 0;
	int nu = 0;
	int i;

	gretl_error_clear();

	nd += seek_file_collections(gretl_home(), DATA_SEARCH, &derr[0]);
	ns += seek_file_collections(gretl_home(), SCRIPT_SEARCH, &serr[0]);
#ifdef OS_OSX
	nd += seek_file_collections(gretl_app_support_dir(), DATA_SEARCH, &derr[1]);
	ns += seek_file_collections(gretl_app_support_dir(), SCRIPT_SEARCH, &serr[1]);
#else
	nd += seek_file_collections(gretl_dotdir(), DATA_SEARCH, &derr[1]);
	ns += seek_file_collections(gretl_dotdir(), SCRIPT_SEARCH, &serr[1]);
#endif

	nu += seek_file_collections(gretl_workdir(), USER_SEARCH, &uerr[0]);
	wd = maybe_get_default_workdir();
	if (wd != NULL) {
	    nu += seek_file_collections(wd, USER_SEARCH, &uerr[1]);
	    nd += seek_file_collections(wd, DATA_SEARCH, &derr[2]);
	    ns += seek_file_collections(wd, SCRIPT_SEARCH, &serr[2]);
	}

	for (i=0; i<3; i++) {
	    /* help to diagnose any errors? */
	    if (derr[i]) {
		fprintf(stderr, "data seek %d gave error %d\n", i+1, derr[i]);
	    }
	    if (serr[i]) {
		fprintf(stderr, "script seek %d gave error %d\n", i+1, serr[i]);
	    }
	    if (i < 2 && uerr[i]) {
		fprintf(stderr, "user file seek %d gave error %d\n", i+1, uerr[i]);
	    }
	}
	if (role == TEXTBOOK_DATA) {
	    if (nd == 0) {
		gretl_errmsg_ensure("file_collections: no data files found");
		errbox(gretl_errmsg_get());
		err = 1;
	    }
	} else if (role == PS_FILES) {
	    if (ns == 0) {
		gretl_errmsg_ensure("file_collections: no script files found");
		errbox(gretl_errmsg_get());
		err = 1;
	    }
	}
	if (!err) {
	    if (nd > 0) {
		sort_files_stack(TEXTBOOK_DATA);
	    }
	    if (ns > 0) {
		sort_files_stack(PS_FILES);
	    }
	}
	collections_built = TRUE;
    }

#if COLL_DEBUG
    print_collections(TEXTBOOK_DATA);
    print_collections(PS_FILES);
#endif

    return err;
}

char *strip_extension (char *s)
{
    char *p = strstr(s, ".tar.gz");

    if (p != NULL) {
	char *q = strstr(s, "_data");

	if (q != NULL) {
	    *q = '\0';
	} else {
	    *p = '\0';
	}
    } else {
	p = strrchr(s, '.');

	if (p != NULL &&
	    (!strcmp(p, ".gdt") || !strcmp(p, ".inp") ||
	     !strcmp(p, ".bin") || !strcmp(p, ".gfn") ||
	     !strcmp(p, ".bn7") || !strcmp(p, ".zip"))) {
	    *p = '\0';
	}
    }

    return s;
}

static int validate_desc_strings (const char *s1,
				  const char *s2,
				  const char *s3)
{
    int err = 0;

    if (!g_utf8_validate(s1, -1, NULL)) {
	err = E_DATA;
    } else if (!g_utf8_validate(s2, -1, NULL)) {
	err = E_DATA;
    } else if (s3 != NULL && !g_utf8_validate(s3, -1, NULL)) {
	err = E_DATA;
    }

    return err;
}

static int read_file_descriptions (windata_t *win, gpointer p)
{
    file_collection *collection = (file_collection *) p;
    GtkListStore *store;
    GtkTreeSelection *selection;
    GtkTreeIter iter;
    char line[MAXLEN];
    char *index;
    FILE *fp;
    int datacols = 0;
    int err = 0;

    index = full_path(collection->path, collection->descfile);

    fp = gretl_fopen(index, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(win->listbox)));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    while (fgets(line, sizeof line, fp) && !err) {
	char fname[24], descrip[80], data[64];
	int nf;

	if (*line == '#') continue;

	nf = sscanf(line, " \"%23[^\"]\",\"%79[^\"]\",\"%63[^\"]\"",
		    fname, descrip, data);

	if (nf == 3) {
	    err = validate_desc_strings(fname, descrip, data);
	    if (!err) {
		datacols = 3;
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter,
				   0, strip_extension(fname),
				   1, descrip,
				   2, data, -1);
	    }
	} else if (nf == 2) {
	    err = validate_desc_strings(fname, descrip, NULL);
	    if (!err) {
		datacols = 2;
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter,
				   0, strip_extension(fname),
				   1, descrip, -1);
	    }
	} else {
	    ; /* ?? */
	}
    }

    fclose(fp);

    if (!err && datacols == 2) {
	/* these windows can have either 2 or 3 columns */
	GtkTreeViewColumn *col;

	col = gtk_tree_view_get_column(GTK_TREE_VIEW(win->listbox), 2);
	if (col != NULL) {
	    gtk_tree_view_column_set_visible(col, FALSE);
	}
    }

    if (!err) {
	/* select the first row */
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(win->listbox));
	gtk_tree_selection_select_iter(selection, &iter);
    }

    return err;
}

static void show_datafile_info (GtkWidget *w, gpointer data)
{
    char fullname[MAXLEN];
    windata_t *vwin = (windata_t *) data;
    file_collection *collection;
    char *descrip;
    gchar *filename;
    int err = 0;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &filename);
    collection = g_object_get_data(G_OBJECT(vwin->listbox), "collection");
    gretl_build_path(fullname, collection->path, filename, NULL);
    strcat(fullname, ".gdt");
    g_free(filename);

#if 0
    fprintf(stderr, "info: active=%d, fullname='%s'\n", vwin->active_var,
	    fullname);
    fprintf(stderr, "collection path='%s'\n", collection->path);
#endif

    descrip = gretl_get_gdt_description(fullname, &err);

    if (err) {
	gui_errmsg(err);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", _("data info"));
	PRN *prn;

	prn = gretl_print_new_with_buffer(descrip);
	view_buffer(prn, 80, 320, title, INFO, NULL);
	g_free(title);
    }
}

void browser_open_data (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    file_collection *collection;
    char tmp[MAXLEN];
    gchar *filename;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &filename);
    collection = g_object_get_data(G_OBJECT(vwin->listbox), "collection");
    gretl_build_path(tmp, collection->path, filename, NULL);
    strcat(tmp, ".gdt");
    set_tryfile(tmp);
    g_free(filename);

    set_datapage(collection->title);

    verify_open_data(vwin, OPEN_DATA, FALSE);
}

void browser_open_ps (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    file_collection *collection;
    gchar *filename;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &filename);
    collection = g_object_get_data(G_OBJECT(vwin->listbox), "collection");
    gretl_build_path(scriptfile, collection->path, filename, NULL);
    strcat(scriptfile, ".inp");
    g_free(filename);

    /* close the calling window */
    gtk_widget_destroy(GTK_WIDGET(vwin->main));

    set_scriptpage(collection->title);

    view_script(scriptfile, 0, VIEW_SCRIPT);
}

static gint enter_opens_file (GtkWidget *w, GdkEventKey *key,
			      windata_t *vwin)
{
    if (key->keyval == GDK_Return || key->keyval == GDK_o) {
	if (vwin->role == TEXTBOOK_DATA) {
	    browser_open_data(w, vwin);
	} else if (vwin->role == PS_FILES) {
	    browser_open_ps(w, vwin);
	}
	return TRUE;
    } else {
	return FALSE;
    }
}

static void browser_delete_current_row (windata_t *vwin)
{
    GtkTreeModel *mod;
    GtkTreeIter iter;
    int i;

    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
    if (!gtk_tree_model_get_iter_first(mod, &iter)) {
	return;
    }

    for (i=0; ; i++) {
	if (i == vwin->active_var) {
	    gtk_list_store_remove(GTK_LIST_STORE(mod), &iter);
	    break;
	} else if (!gtk_tree_model_iter_next(mod, &iter)) {
	    break;
	}
    }
}

static void browser_delete_row_by_content (windata_t *vwin,
					   int colnum1,
					   const char *test1,
					   int colnum2,
					   const char *test2)
{
    GtkTreeModel *mod;
    GtkTreeIter iter;
    gchar *content;
    int done = 0;

    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
    if (!gtk_tree_model_get_iter_first(mod, &iter)) {
	return;
    }

    while (1) {
	gtk_tree_model_get(mod, &iter, colnum1, &content, -1);
	if (content != NULL && !strcmp(content, test1)) {
	    if (test2 == NULL) {
		done = 1;
	    } else {
		g_free(content);
		gtk_tree_model_get(mod, &iter, colnum2, &content, -1);
		if (content != NULL && !strcmp(content, test2)) {
		    done = 1;
		}
	    }
	    if (done) {
		gtk_list_store_remove(GTK_LIST_STORE(mod), &iter);
	    }
	}
	g_free(content);
	if (done || !gtk_tree_model_iter_next(mod, &iter)) {
	    break;
	}
    }
}

static int gui_delete_fn_pkg (const char *pkgname, const char *fname,
			      windata_t *vwin)
{
    gchar *msg = NULL;
    int delete_ok = 0;
    int delete = 0;
    int loaded = 0;
    int unload = 0;
    int resp, err = 0;

    if (package_being_edited(pkgname, NULL)) {
	warnbox_printf(_("%s: please close this object's window first"),
		       pkgname);
	return 0;
    }

    /* see if the package is loaded in memory */
    if (get_function_package_by_name(pkgname) != NULL) {
	loaded = 1;
    }

    /* see if the user is able to delete the package */
    if (!is_gretl_addon(pkgname) && gretl_write_access((char *) fname) == 0) {
	delete_ok = 1;
    }

    if (!loaded && !delete_ok) {
	infobox_printf(_("Package %s is not loaded, and you do "
			 "not have permission to delete it."),
		       pkgname);
	return 0;
    }

    if (loaded && delete_ok) {
	const char *opts[] = {
	    N_("Unload member functions only"),
	    N_("Unload and delete package file"),
	};

	msg = g_strdup_printf(_("Function package %s"), fname);
	resp = radio_dialog(NULL, msg, opts, 2, 0, 0, vwin_toplevel(vwin));
	g_free(msg);
	if (resp < 0) {
	    /* canceled */
	    return 0;
	} else if (resp == 1) {
	    delete = 1;
	} else {
	    unload = 1;
	}
    } else if (loaded) {
	msg = g_strdup_printf(_("Unload package %s?"), pkgname);
	resp = yes_no_dialog(NULL, msg, vwin_toplevel(vwin));
	g_free(msg);
	if (resp == GRETL_NO) {
	    return 0;
	} else {
	    unload = 1;
	}
    } else if (delete_ok) {
	msg = g_strdup_printf(_("Really delete %s?"), pkgname);
	resp = yes_no_dialog(NULL, msg, vwin_toplevel(vwin));
	g_free(msg);
	if (resp == GRETL_NO) {
	    return 0;
	} else {
	    delete = 1;
	}
    }

    if (unload && !delete) {
	/* just unload the package from memory */
	function_package_unload_full_by_filename(fname);
    } else if (delete) {
	/* remove entry from registry, if present */
	gui_function_pkg_unregister(pkgname);
	/* unload the package from memory */
	if (loaded) {
	    function_package_unload_full_by_filename(fname);
	}
	/* trash the package file(s) */
	err = delete_function_package(fname);
	if (err) {
	    gui_errmsg(err);
	} else {
	    /* remove package from GUI listing */
	    browser_delete_current_row(vwin);
	}
    }

    return err;
}

static int get_info_width (const char *buf)
{
    char line[1024];
    int n, width = 66;

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	n = strlen(line);
	if (n > width && n <= HELP_WIDTH) {
	    width = n;
	}
    }

    bufgets_finalize(buf);

    return width + 1;
}

windata_t *display_function_package_data (const char *pkgname,
					  const char *path,
					  int role)
{
    windata_t *vwin = NULL;
    PRN *prn = NULL;
    int err = 0;

    if (bufopen(&prn)) {
	return NULL;
    }

    if (role == VIEW_PKG_INFO) {
	err = print_function_package_info(path, 1, prn);
    } else if (role == VIEW_PKG_SAMPLE) {
	err = print_function_package_sample(path, tabwidth, prn);
    } else {
	err = print_function_package_code(path, tabwidth, prn);
    }

    if (err) {
	gretl_print_destroy(prn);
	gui_errmsg(err);
    } else {
	gchar *title;

	if (role == VIEW_PKG_SAMPLE) {
	    title = g_strdup_printf("gretl: %s sample", pkgname);
	} else {
	    title = g_strdup_printf("gretl: %s", pkgname);
	}

	if (role == VIEW_PKG_INFO) {
	    char *buf = gretl_print_steal_buffer(prn);
	    int width = get_info_width(buf);

	    vwin = view_formatted_text_buffer(title, buf, width, 350, role);
	    free(buf);
	    gretl_print_destroy(prn);
	} else {
	    vwin = view_buffer(prn, 78, 350, title, role, NULL);
	}
	strcpy(vwin->fname, path);
	if (strstr(path, "dltmp")) {
	    set_window_delete_filename(vwin);
	}
	g_free(title);
    }

    return vwin;
}

/* on adding files to gfn view window: ensure that the GTK
   selection stays in sync with the "active_var" ID
*/

static void fix_selected_row (GtkTreeModel *model,
			      GtkTreePath *path,
			      GtkTreeIter *iter,
			      gpointer data)
{
    gint idx = gtk_tree_path_get_indices(path)[0];
    windata_t *vwin = data;

    vwin->active_var = idx;
}

/* callback from the file selector where the user has chosen
   a directory at which to point the gfn browser
*/

void set_alternate_gfn_dir (windata_t *vwin, char *path)
{
    GDir *dir;
    int replace = 1;
    int nfn = 0;

#if GFN_DEBUG
    fprintf(stderr, "set_alternate_gfn_dir: '%s'\n", path);
#endif

    dir = gretl_opendir(path);
    if (dir == NULL) {
	/* should never happen, but... */
	return;
    }

    /* first pass: just count gfn files in @path */
    read_fn_files_in_dir(dir, path, NULL, NULL, &nfn);

    if (nfn == 0) {
	warnbox(_("No function files were found"));
	replace = 0;
    } else {
	/* give the user a chance to back out */
	gchar *msg;

	msg = g_strdup_printf(_("Found %d function file(s).\n"
				"Replace the current listing?"),
			      nfn);
	if (yes_no_dialog("gretl", msg, vwin->main) != GRETL_YES) {
	    replace = 0;
	}
	g_free(msg);
    }

    if (replace) {
	/* OK: now rewind, clear, and repopulate the listing using
	   the @path selected by the user
	*/
	GtkListStore *store;
	GtkTreeSelection *sel;
	GtkTreeIter iter;
	gulong sigid;

	store = GTK_LIST_STORE(gtk_tree_view_get_model
			       (GTK_TREE_VIEW(vwin->listbox)));
	g_dir_rewind(dir);
	sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
	sigid = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sel), "menu_check"));
	g_signal_handler_block(sel, sigid);
	gtk_list_store_clear(store);
	g_signal_handler_unblock(sel, sigid);
	nfn = 0;
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	read_fn_files_in_dir(dir, path, store, &iter, &nfn);

	if (nfn > 0) {
	    gtk_tree_selection_selected_foreach(sel, fix_selected_row, vwin);
	    widget_set_int(vwin->listbox, "altdir", 1);
	    presort_treelist(vwin, NULL);
	    listbox_select_first(vwin);
	} else {
	    /* can't happen? */
	    warnbox(_("No function files were found"));
	}
    }

    g_dir_close(dir);
}

gchar *gfn_browser_get_alt_path (void)
{
    windata_t *vwin = get_browser_for_role(FUNC_FILES, NULL);
    gchar *path = NULL;

    if (vwin != NULL && widget_get_int(vwin->listbox, "altdir")) {
	tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 0,
			     GFN_DIRNAME_COL, &path);
    }

    return path;
}

static void query_remove_gfn_from_registry (const char *pkgname,
					    windata_t *vwin)
{
    gchar *msg;
    int resp;

    msg = g_strdup_printf(_("Really remove %s from menu?"), pkgname);
    resp = yes_no_dialog(NULL, msg, vwin->main);

    if (resp == GRETL_YES) {
	gui_function_pkg_unregister(pkgname);
    }

    g_free(msg);
}

static void browser_functions_handler (windata_t *vwin, int task)
{
    char path[FILENAME_MAX];
    gchar *pkgname = NULL;
    gchar *dir = NULL;
    int dircol = 0;

    if (vwin->role == DBNOMICS_TOP && task == VIEW_PKG_DOC) {
	/* special case: not coming from gfn window */
	gchar *docpath = g_build_filename(gretl_home(), "functions",
					  "dbnomics", "dbnomics.pdf",
					  NULL);
	gretl_show_pdf(docpath, NULL);
	g_free(docpath);
	return;
    }

    if (vwin->role == FUNC_FILES) {
	dircol = GFN_DIRNAME_COL;
    } else if (vwin->role != REMOTE_FUNC_FILES &&
	       vwin->role != PKG_REGISTRY) {
	dummy_call();
	return;
    }

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &pkgname);

    if (dircol != 0) {
	tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			     dircol, &dir);
	if (task == VIEW_PKG_RESOURCES) {
	    gretl_build_path(path, dir, "examples", NULL);
	} else if (task == VIEW_PKG_DOC) {
	    gretl_build_path(path, dir, pkgname, NULL);
	    strcat(path, ".pdf");
	} else {
	    /* VIEW_FN_PKG_INFO */
	    gretl_build_path(path, dir, pkgname, NULL);
	    strcat(path, ".gfn");
	}
    } else {
	strcpy(path, pkgname);
    }

    if (vwin->role == FUNC_FILES && task != DELETE_FN_PKG) {
	/* try to ensure we don't get a stale version */
	char test[FILENAME_MAX];
	const char *loaded = NULL;
	gchar *version = NULL;

	gretl_build_path(test, dir, pkgname, NULL);
	strcat(test, ".gfn");
	if (function_package_is_loaded(test, &loaded)) {
	    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
				 1, &version);
	    if (loaded != NULL && version != NULL &&
		strcmp(version, loaded)) {
		function_package_unload_by_filename(test);
	    }
	}
	g_free(version);
    }

#if GFN_DEBUG
    fprintf(stderr, "browser_functions_handler: active=%d, pkgname='%s'\n"
	    "path='%s'\n", vwin->active_var, pkgname, path);
#endif

    if (task == DELETE_FN_PKG) {
	gui_delete_fn_pkg(pkgname, path, vwin);
    } else if (task == VIEW_FN_PKG_INFO) {
	display_function_package_data(pkgname, path, VIEW_PKG_INFO);
    } else if (task == VIEW_FN_PKG_SAMPLE) {
	display_function_package_data(pkgname, path, VIEW_PKG_SAMPLE);
    } else if (task == VIEW_FN_PKG_CODE) {
	display_function_package_data(pkgname, path, VIEW_PKG_CODE);
    } else if (task == MENU_ADD_FN_PKG) {
	gui_function_pkg_query_register(path, vwin->main);
    } else if (task == MENU_REMOVE_FN_PKG) {
	query_remove_gfn_from_registry(pkgname, vwin);
    } else if (task == VIEW_PKG_RESOURCES) {
	file_selector_with_startdir(OPEN_ANY, path, vwin_toplevel(vwin));
    } else if (task == VIEW_PKG_DOC) {
	gretl_show_pdf(path, NULL);
    } else if (task == CALL_FN_PKG) {
	/* note: this is the double-click default for the local
	   function package browser; can also be invoked via
	   context menu and toolbar button
	*/
	open_function_package(pkgname, path, vwin);
    }

    g_free(pkgname);
    g_free(dir);
}

static void show_addon_info (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    int v = vwin->active_var;
    gchar *pkgname = NULL;
    gchar *pkgpath = NULL;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), v, 0, &pkgname);
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), v, 3, &pkgpath);

    if (pkgname != NULL && pkgpath != NULL) {
	display_function_package_data(pkgname, pkgpath, VIEW_PKG_INFO);
	g_free(pkgname);
	g_free(pkgpath);
    } else {
	errbox("Couldn't get info for addon");
    }
}

/* this function is public because it's called from
   doubleclick_action() in callbacks.c
*/

void browser_call_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, CALL_FN_PKG);
}

static void show_function_info (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, VIEW_FN_PKG_INFO);
}

static void show_function_code (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, VIEW_FN_PKG_CODE);
}

static void show_function_sample (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, VIEW_FN_PKG_SAMPLE);
}

static void show_package_resources (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, VIEW_PKG_RESOURCES);
}

static void show_package_doc (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, VIEW_PKG_DOC);
}

static void browser_del_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, DELETE_FN_PKG);
}

static void add_func_to_menu (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, MENU_ADD_FN_PKG);
}

static void gfn_registry_remove (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, MENU_REMOVE_FN_PKG);
}

windata_t *get_local_viewer (int remote_role)
{
    windata_t *vwin = NULL;

    if (remote_role == REMOTE_DB) {
	vwin = get_browser_for_role(NATIVE_DB, NULL);
    } else if (remote_role == REMOTE_FUNC_FILES) {
	vwin = get_browser_for_role(FUNC_FILES, NULL);
    }

    return vwin;
}

void start_new_function_package (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    functions_selection_wrapper(vwin_toplevel(vwin));
}

static void build_datafiles_popup (windata_t *vwin)
{
    if (vwin->popup == NULL) {
	vwin->popup = gtk_menu_new();
	add_popup_item(_("Info"), vwin->popup,
		       G_CALLBACK(show_datafile_info),
		       vwin);
	add_popup_item(_("Open"), vwin->popup,
		       G_CALLBACK(browser_open_data),
		       vwin);
    }
}

static int get_menu_add_ok (windata_t *vwin)
{
    gchar *pkgname = NULL;
    gchar *dirname = NULL;
    int ret = 0;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &pkgname);
    if (pkgname != NULL && !strcmp(pkgname, "ridge")) {
	return 0;
    }
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 GFN_DIRNAME_COL, &dirname);

    if (pkgname != NULL && dirname != NULL) {
	char path[FILENAME_MAX];

	gretl_build_path(path, dirname, pkgname, NULL);
	strcat(path, ".gfn");
	ret = package_is_available_for_menu(pkgname, path);
    }

    g_free(pkgname);
    g_free(dirname);

    return ret;
}

static void check_extra_buttons_state (GtkTreeSelection *sel, windata_t *vwin)
{
    GtkWidget *button;
    gint flags = 0;

    button = g_object_get_data(G_OBJECT(vwin->mbar), "add-button");
    if (button != NULL) {
	gtk_widget_set_sensitive(button, get_menu_add_ok(vwin));
    }

    /* Get flags from last, hidden int column, to determine
       whether we can offer links to a package's "examples"
       directory and/or its documentation in PDF format.
    */
    tree_view_get_int(GTK_TREE_VIEW(vwin->listbox),
		      vwin->active_var, GFN_FLAGS_COL, &flags);

    button = g_object_get_data(G_OBJECT(vwin->mbar), "res-button");
    if (button != NULL) {
	gtk_widget_set_sensitive(button, flags & PKG_ATTR_RES);
    }

    button = g_object_get_data(G_OBJECT(vwin->mbar), "doc-button");
    if (button != NULL) {
	gtk_widget_set_sensitive(button, flags & PKG_ATTR_DOC);
    }
}

static void connect_menu_adjust_signal (windata_t *vwin)
{
    GtkTreeSelection *sel;
    gulong id;

    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
    id = g_signal_connect(G_OBJECT(sel), "changed",
			  G_CALLBACK(check_extra_buttons_state),
			  vwin);
    g_object_set_data(G_OBJECT(sel), "menu_check", GINT_TO_POINTER(id));
}

static void build_funcfiles_popup (windata_t *vwin)
{
    vwin->popup = gtk_menu_new();

    if (vwin->role == FUNC_FILES) {
	/* local function files: full menu */
	int add_ok = 0;
	int res_ok = 0;
	int doc_ok = 0;
	GtkWidget *b;

	b = g_object_get_data(G_OBJECT(vwin->mbar), "add-button");
	if (b != NULL && gtk_widget_is_sensitive(b)) {
	    add_ok = 1;
	}

	b = g_object_get_data(G_OBJECT(vwin->mbar), "res-button");
	if (b != NULL && gtk_widget_is_sensitive(b)) {
	    res_ok = 1;
	}

	b = g_object_get_data(G_OBJECT(vwin->mbar), "doc-button");
	if (b != NULL && gtk_widget_is_sensitive(b)) {
	    doc_ok = 1;
	}

	add_popup_item(_("Info"), vwin->popup,
		       G_CALLBACK(show_function_info),
		       vwin);
	add_popup_item(_("Sample script"), vwin->popup,
		       G_CALLBACK(show_function_sample),
		       vwin);
	add_popup_item(_("View code"), vwin->popup,
		       G_CALLBACK(show_function_code),
		       vwin);
	add_popup_item(_("Execute"), vwin->popup,
		       G_CALLBACK(browser_call_func),
		       vwin);
	if (res_ok) {
	    add_popup_item(_("Resources..."), vwin->popup,
			   G_CALLBACK(show_package_resources),
			   vwin);
	}
	if (add_ok) {
	    add_popup_item(_("Add to menu"), vwin->popup,
			   G_CALLBACK(add_func_to_menu),
			   vwin);
	}
	if (doc_ok) {
	    add_popup_item(_("Help"), vwin->popup,
			   G_CALLBACK(show_package_doc),
			   vwin);
	}
	add_popup_item(_("Unload/delete..."), vwin->popup,
		       G_CALLBACK(browser_del_func),
		       vwin);
    } else if (vwin->role == REMOTE_FUNC_FILES) {
	/* files on server: limited menu */
	add_popup_item(_("Info"), vwin->popup,
		       G_CALLBACK(pkg_info_from_server),
		       vwin);
	add_popup_item(_("Install"), vwin->popup,
		       G_CALLBACK(install_file_from_server),
		       vwin);
    } else if (vwin->role == ADDONS_FILES) {
	add_popup_item(_("Info"), vwin->popup,
		       G_CALLBACK(show_addon_info),
		       vwin);
    } else if (vwin->role == PKG_REGISTRY) {
	add_popup_item(_("Remove"), vwin->popup,
		       G_CALLBACK(gfn_registry_remove),
		       vwin);
    }
}

static gboolean
funcfiles_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer data)
{
    if (right_click(event)) {
	windata_t *vwin = (windata_t *) data;

	if (vwin->popup != NULL) {
	    gtk_widget_destroy(vwin->popup);
	    vwin->popup = NULL;
	}

	build_funcfiles_popup(vwin);

	if (vwin->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(vwin->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    g_signal_connect(G_OBJECT(vwin->popup), "destroy",
			     G_CALLBACK(gtk_widget_destroyed),
			     &vwin->popup);
	}

	return TRUE;
    }

    return FALSE;
}

static void build_db_popup (windata_t *vwin)
{
    if (vwin->popup == NULL) {
	vwin->popup = gtk_menu_new();
	if (vwin->role == NATIVE_DB) {
	    add_popup_item(_("List series"), vwin->popup,
			   G_CALLBACK(open_db_index),
			   vwin);
	    add_popup_item(_("Find..."), vwin->popup,
			   G_CALLBACK(listbox_find),
			   vwin);
	} else {
	    add_popup_item(_("List series"), vwin->popup,
			   G_CALLBACK(open_remote_db_index),
			   vwin);
	    add_popup_item(_("Install"), vwin->popup,
			   G_CALLBACK(install_file_from_server),
			   vwin);
	    add_popup_item(_("Find..."), vwin->popup,
			   G_CALLBACK(listbox_find),
			   vwin);
	}
    }
}

static void build_data_pkg_popup (windata_t *vwin)
{
    if (vwin->popup == NULL) {
	vwin->popup = gtk_menu_new();
	add_popup_item(_("Install"), vwin->popup,
		       G_CALLBACK(install_file_from_server),
		       vwin);
    }
}

static void show_server_dbs (GtkWidget *w, gpointer p)
{
    display_files(REMOTE_DB, NULL);
}

static void show_local_dbs (GtkWidget *w, gpointer p)
{
    display_files(NATIVE_DB, NULL);
}

static void show_server_funcs (GtkWidget *w, gpointer p)
{
    display_files(REMOTE_FUNC_FILES, NULL);
}

static void show_server_data_pkgs (GtkWidget *w, gpointer p)
{
    display_files(REMOTE_DATA_PKGS, NULL);
}

static void show_local_funcs (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, NULL);
}

static void show_gfn_registry (GtkWidget *w, windata_t *vwin)
{
    display_files(PKG_REGISTRY, NULL);
}

/* Respond when the user has clicked the Directory button
   in the function package browser. What exactly we do
   here depends on whether the browser is currently in
   its default mode (viewing installed packages) or if
   it is redirected -- which is flagged by a non-zero
   value for "altdir" on the browser's listbox.
*/

static void alt_funcs_dir (GtkWidget *w, windata_t *vwin)
{
    if (widget_get_int(vwin->listbox, "altdir")) {
	const char *opts[] = {
	    N_("Choose another directory"),
	    N_("Revert to installed packages")
	};
	int resp;

	resp = radio_dialog(NULL, NULL, opts, 2, 0, 0, vwin->main);

	if (resp == GRETL_CANCEL) {
	    return;
	} else if (resp == 1) {
	    /* revert to installed gfns */
	    widget_set_int(vwin->listbox, "altdir", 0);
	    populate_gfn_list(vwin);
	    listbox_select_first(vwin);
	    return;
	}
    }

    /* If not canceled or reverted to default, put up a
       dialog to let the user select a directory: the
       callback from that is set_alternate_gfn_dir().
    */

    file_selector_with_parent(SET_FDIR, FSEL_DATA_VWIN, vwin,
                              vwin->main);
}

static void alt_db_dir (GtkWidget *w, windata_t *vwin)
{
    file_selector_with_parent(SET_DBDIR, FSEL_DATA_VWIN, vwin,
                              vwin->main);
}

static GtkWidget *get_dbn_menu (windata_t *vwin)
{
    GtkWidget *menu = gtk_menu_new();
    GtkAction *action;
    GtkWidget *item;

    action = gtk_action_new("DBNbrowse", _("Browse..."), NULL, NULL);
    g_signal_connect(G_OBJECT(action), "activate",
		     G_CALLBACK(show_files), vwin);
    item = gtk_action_create_menu_item(action);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);

    action = gtk_action_new("DBNseries", _("Specific series..."), NULL, NULL);
    g_signal_connect(G_OBJECT(action), "activate",
		     G_CALLBACK(dbnomics_specific_series), vwin);
    item = gtk_action_create_menu_item(action);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);

    /* don't leak: record pointer to menu so it can
       be destroyed when the window is closed */
    vwin_record_toolbar_popup(vwin, menu);

    return menu;
}

enum {
    BTN_INFO = 1,
    BTN_CODE,
    BTN_INDX,
    BTN_INST,
    BTN_EXEC,
    BTN_ADD,
    BTN_DEL,
    BTN_WWW,
    BTN_HOME,
    BTN_OPEN,
    BTN_DIR,
    BTN_RES,
    BTN_DOC,
    BTN_REG,
    BTN_DBN,
    BTN_PROP
};

static GretlToolItem files_items[] = {
    { N_("Open"),           GTK_STOCK_OPEN, NULL, BTN_OPEN },
    { N_("Select directory"), GTK_STOCK_DIRECTORY, NULL, BTN_DIR },
    { N_("Info"),           GTK_STOCK_INFO,       NULL, BTN_INFO },
    { N_("Details"),        GTK_STOCK_PROPERTIES, G_CALLBACK(open_dbnomics_series), BTN_PROP },
    { N_("Sample script"),  GTK_STOCK_JUSTIFY_LEFT, G_CALLBACK(show_function_sample), BTN_CODE },
    { N_("View code"),      GTK_STOCK_PROPERTIES, G_CALLBACK(show_function_code), BTN_CODE },
    { N_("Execute"),        GTK_STOCK_EXECUTE,    G_CALLBACK(browser_call_func),  BTN_EXEC },
    { N_("List series"),    GTK_STOCK_INDEX,      NULL,                           BTN_INDX },
    { N_("Install"),        GTK_STOCK_SAVE,       NULL,                           BTN_INST },
    { N_("Resources..."),   GTK_STOCK_OPEN,       G_CALLBACK(show_package_resources), BTN_RES },
    { N_("Add to menu"),    GTK_STOCK_ADD,        G_CALLBACK(add_func_to_menu),  BTN_ADD },
    { N_("Package registry"), GTK_STOCK_PREFERENCES, G_CALLBACK(show_gfn_registry), BTN_REG },
    { N_("Help"),           GRETL_STOCK_PDF,      G_CALLBACK(show_package_doc),  BTN_DOC },
    { N_("Unload/delete..."), GTK_STOCK_DELETE,   G_CALLBACK(browser_del_func),  BTN_DEL },
    { N_("Look on server"), GTK_STOCK_NETWORK,    NULL,                          BTN_WWW },
    { N_("Local machine"),  GTK_STOCK_HOME,       NULL,                          BTN_HOME },
    { "DB.NOMICS",          GRETL_STOCK_DBN,      NULL,                          BTN_DBN }
};

static GretlToolItem pager_items[] = {
    { N_("First"),    GTK_STOCK_GOTO_FIRST, G_CALLBACK(dbnomics_pager_call), 1 },
    { N_("Previous"), GTK_STOCK_GO_BACK,    G_CALLBACK(dbnomics_pager_call), 2 },
    { N_("Next"),     GTK_STOCK_GO_FORWARD, G_CALLBACK(dbnomics_pager_call), 3 },
    { N_("Last"),     GTK_STOCK_GOTO_LAST,  G_CALLBACK(dbnomics_pager_call), 4 }
};

static int n_files_items = G_N_ELEMENTS(files_items);

#define common_item(f) (f == 0)

#define local_funcs_item(f) (f == BTN_DEL || f == BTN_CODE || \
			     f == BTN_RES || f == BTN_DOC || \
			     f == BTN_REG)

static int files_item_get_callback (GretlToolItem *item, int role)
{
    if (common_item(item->flag)) {
	return 1;
    } else if (item->flag == BTN_DEL) {
	if (role == PKG_REGISTRY) {
	    item->func = G_CALLBACK(gfn_registry_remove);
	    item->tip = N_("Remove from menu");
	    return 1;
	} else if (role == FUNC_FILES) {
	    item->func = G_CALLBACK(browser_del_func);
	    item->tip = N_("Unload/delete...");
	    return 1;
	}
    } else if (local_funcs_item(item->flag)) {
	if (item->flag == BTN_DOC && role == DBNOMICS_TOP) {
	    return 1;
	} else {
	    return (role == FUNC_FILES);
	}
    } else if (item->flag == BTN_INST) {
	if (role == ADDONS_FILES) {
	    return 0;
	} else {
	    item->func = G_CALLBACK(install_file_from_server);
	    return (role == REMOTE_DB ||
		    role == REMOTE_FUNC_FILES ||
		    role == REMOTE_DATA_PKGS);
	}
    } else if (item->flag == BTN_EXEC || item->flag == BTN_ADD) {
	return (role == FUNC_FILES);
    } else if (item->flag == BTN_PROP) {
	return (role == DBNOMICS_SERIES);
    }

    item->func = NULL;

    if (item->flag == BTN_OPEN) {
	/* open: only data files and scripts */
	if (role == TEXTBOOK_DATA) {
	    item->func = G_CALLBACK(browser_open_data);
	} else if (role == PS_FILES) {
	    item->func = G_CALLBACK(browser_open_ps);
	}
    } else if (item->flag == BTN_INFO) {
	if (role == TEXTBOOK_DATA) {
	    item->func = G_CALLBACK(show_datafile_info);
	} else if (role == FUNC_FILES) {
	    item->func = G_CALLBACK(show_function_info);
	} else if (role == REMOTE_FUNC_FILES) {
	    item->func = G_CALLBACK(pkg_info_from_server);
	} else if (role == ADDONS_FILES) {
	    item->func = G_CALLBACK(show_addon_info);
	} else if (role == DBNOMICS_DB) {
	    item->func = G_CALLBACK(show_dbnomics_dimensions);
	}
    } else if (item->flag == BTN_INDX) {
	/* index: databases only */
	item->tip = N_("List series");
	if (role == NATIVE_DB) {
	    item->func = G_CALLBACK(open_db_index);
	} else if (role == REMOTE_DB) {
	    item->func = G_CALLBACK(open_remote_db_index);
	} else if (role == DBNOMICS_TOP) {
	    item->func = G_CALLBACK(open_dbnomics_provider);
	    item->tip = N_("List datasets");
	} else if (role == DBNOMICS_DB) {
	    item->func = G_CALLBACK(open_dbnomics_dataset);
	}
    } else if (item->flag == BTN_WWW) {
	if (role == FUNC_FILES) {
	    item->func = G_CALLBACK(show_server_funcs);
	} else if (role == NATIVE_DB) {
	    item->func = G_CALLBACK(show_server_dbs);
	} else if (role == TEXTBOOK_DATA) {
	    item->func = G_CALLBACK(show_server_data_pkgs);
	}
    } else if (item->flag == BTN_HOME) {
	/* home: show only for on-server items */
	if (role == REMOTE_FUNC_FILES) {
	    item->func = G_CALLBACK(show_local_funcs);
	} else if (role == REMOTE_DB) {
	    item->func = G_CALLBACK(show_local_dbs);
	}
    } else if (item->flag == BTN_DIR) {
	if (role == FUNC_FILES) {
	    item->func = G_CALLBACK(alt_funcs_dir);
	} else if (role == NATIVE_DB) {
	    item->func = G_CALLBACK(alt_db_dir);
	}
    } else if (item->flag == BTN_DBN) {
	if (role == NATIVE_DB) {
	    item->func = G_CALLBACK(dummy_call);
	}
    }

    return (item->func != NULL);
}

static void filter_remote_funcs (GtkButton *b, windata_t *vwin)
{
    GtkWidget *combo;
    int filter, old_filter;
    gchar *s;

    combo = g_object_get_data(G_OBJECT(vwin->main), "filter-combo");
    old_filter = widget_get_int(combo, "filter");

    if (gtk_combo_box_get_active(GTK_COMBO_BOX(combo)) == 0) {
	filter = 0;
    } else {
	s = combo_box_get_active_text(combo);
	filter = atoi(s + 1); /* "C<number>" */
	g_free(s);
    }

    if (filter != old_filter) {
	populate_remote_func_list(vwin, filter);
	widget_set_int(combo, "filter", filter);
	listbox_select_first(vwin);
    }
}

static void maybe_add_gfn_filter (windata_t *vwin,
				  GtkWidget *hbox)
{
    char *getbuf = NULL;
    int err;

    err = list_remote_function_categories(&getbuf, OPT_NONE);
    if (!err && (getbuf == NULL || *getbuf != 'C')) {
	free(getbuf);
	err = 1;
    }

#if 0
    fprintf(stderr, "getbuf: '%s'\n", getbuf);
#endif

    if (!err) {
	GtkWidget *combo, *button;
	char line[128];
	gchar *label;

	bufgets_init(getbuf);

	combo = gtk_combo_box_text_new();
	g_object_set_data(G_OBJECT(vwin->main), "filter-combo", combo);
	widget_set_int(combo, "filter", 0);

	combo_box_append_text(combo, _("All packages"));
	while (bufgets(line, sizeof line, getbuf)) {
	    combo_box_append_text(combo, tailstrip(line));
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);

	bufgets_finalize(getbuf);
	free(getbuf);

	gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
	label = g_strdup_printf(" %s ", _("filter"));
	button = gtk_button_new_with_label(label);
	g_free(label);
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(filter_remote_funcs), vwin);

    }
}

static void add_dbnomics_menu_button (GretlToolItem *item,
				      windata_t *vwin)
{
    GtkWidget *menu = get_dbn_menu(vwin);

    vwin_toolbar_insert(item, NULL, menu, vwin, -1);
}

static void make_files_toolbar (windata_t *vwin)
{
    GtkWidget *hbox, *button;
    GretlToolItem *item;
    int i;

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gretl_toolbar_new(NULL);

    for (i=0; i<n_files_items; i++) {
	item = &files_items[i];
	if (item->flag == BTN_DBN) {
	    if (vwin->role == NATIVE_DB) {
		add_dbnomics_menu_button(item, vwin);
	    }
	    continue;
	}
	if (files_item_get_callback(item, vwin->role)) {
	    button = gretl_toolbar_insert(vwin->mbar, item, item->func, vwin, -1);
	    if (item->flag == BTN_ADD) {
		g_object_set_data(G_OBJECT(vwin->mbar), "add-button", button);
		gtk_widget_set_sensitive(button, FALSE);
	    } else if (item->flag == BTN_RES) {
		g_object_set_data(G_OBJECT(vwin->mbar), "res-button", button);
		gtk_widget_set_sensitive(button, FALSE);
	    } else if (item->flag == BTN_DOC && vwin->role != DBNOMICS_TOP) {
		/* condition on availablity of PDF documentation */
		g_object_set_data(G_OBJECT(vwin->mbar), "doc-button", button);
		gtk_widget_set_sensitive(button, FALSE);
	    }
	}
    }

    if (DBNOMICS_ACTION(vwin->role)) {
	const char *ids[] = {
	    "first-button", "prev-button", "next-button", "last-button"
	};
	for (i=0; i<4; i++) {
	    item = &pager_items[i];
	    button = gretl_toolbar_insert(vwin->mbar, item, item->func, vwin, -1);
	    g_object_set_data(G_OBJECT(vwin->mbar), ids[i], button);
	    widget_set_int(button, "action", item->flag);
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);

    if (vwin->role == REMOTE_FUNC_FILES) {
	maybe_add_gfn_filter(vwin, hbox);
    }

    vwin_add_winlist(vwin);
    if (vwin->role != ADDONS_FILES) {
	/* there aren't enough addons to warrant this */
	vwin_add_finder(vwin);
    }
    gtk_widget_show_all(hbox);
}

static gchar *files_title (int code)
{
    static char hname[48];
    gchar *ret = NULL;

    if (*hname == '\0') {
	const gchar *s = g_get_host_name();

	if (s != NULL && strlen(s) < 48) {
	    strcpy(hname, s);
	} else {
	    strcpy(hname, _("local machine"));
	}
    }

    if (code == NATIVE_DB) {
	ret = g_strdup_printf(_("gretl: databases on %s"), hname);
    } else {
	ret = g_strdup_printf(_("gretl: function packages on %s"), hname);
    }

    return ret;
}

/* handle drag of pointer from remote database window
   to local one, or remote function package window to
   local one
*/

static void
remote_window_handle_drag (GtkWidget *widget,
			   GdkDragContext *context,
			   gint x,
			   gint y,
			   GtkSelectionData *data,
			   guint info,
			   guint time,
			   gpointer p)
{
    if (data != NULL &&
	(info == GRETL_REMOTE_DB_PTR ||
	 info == GRETL_REMOTE_FNPKG_PTR)) {
	drag_file_from_server(info);
    }
}

static void set_up_viewer_drag_target (windata_t *vwin)
{
    GCallback callback;
    int i;

    if (vwin->role == NATIVE_DB) {
	i = GRETL_REMOTE_DB_PTR;
	callback = G_CALLBACK(remote_window_handle_drag);
    } else if (vwin->role == FUNC_FILES) {
	i = GRETL_REMOTE_FNPKG_PTR;
	callback = G_CALLBACK(remote_window_handle_drag);
    } else {
	return;
    }

    gtk_drag_dest_set(vwin->listbox,
		      GTK_DEST_DEFAULT_ALL,
		      &gretl_drag_targets[i], 1,
		      GDK_ACTION_COPY);

    g_signal_connect(G_OBJECT(vwin->listbox), "drag-data-received",
		     callback, NULL);
}

void listbox_select_first (windata_t *vwin)
{
    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
    GtkTreeModel *model;
    GtkTreeSelection *selection;
    GtkTreeIter iter;

    model = gtk_tree_view_get_model(view);
    gtk_tree_model_get_iter_first(model, &iter);
    selection = gtk_tree_view_get_selection(view);
    gtk_tree_selection_select_iter(selection, &iter);
    gtk_widget_grab_focus(vwin->listbox);
}

static gint catch_browser_key (GtkWidget *w, GdkEventKey *event,
			       windata_t *vwin)
{
    int Ctrl = (event->state & GDK_CONTROL_MASK);
    guint key = event->keyval;

    if (is_control_key(event->keyval) || !Ctrl) {
	return FALSE;
    }
    if (!gdk_keyval_is_upper(event->keyval)) {
	key = gdk_keyval_to_upper(event->keyval);
    }

    if (vwin->role == FUNC_FILES && key == GDK_R) {
	/* "Execute" */
	browser_call_func(NULL, vwin);
	return TRUE;
    } else if (key == GDK_S) {
	/* "Info" */
	if (vwin->role == TEXTBOOK_DATA) {
	    show_datafile_info(NULL, vwin);
	} else if (vwin->role == FUNC_FILES) {
	    show_function_info(NULL, vwin);
	} else if (vwin->role == REMOTE_FUNC_FILES) {
	    pkg_info_from_server(NULL, vwin);
	} else if (vwin->role == ADDONS_FILES) {
	    show_addon_info(NULL, vwin);
	} else if (vwin->role == DBNOMICS_DB) {
	    show_dbnomics_dimensions(NULL, vwin);
	}
	return TRUE;
    }

    return FALSE;
}

static void set_up_browser_keystrokes (windata_t *vwin)
{
    g_signal_connect(G_OBJECT(vwin->main), "key-press-event",
		     G_CALLBACK(catch_browser_key), vwin);
}

#define want_keyvals(r) (r==FUNC_FILES || r==TEXTBOOK_DATA || \
			 r==REMOTE_FUNC_FILES || r==ADDONS_FILES || \
			 r==DBNOMICS_DB)

#define notebook_needed(r) (r==TEXTBOOK_DATA || r==PS_FILES)

void display_files (int role, const gchar *path)
{
    GtkWidget *filebox;
    windata_t *vwin;
    gchar *title = NULL;
    gchar *mypath = NULL;
    int err = 0;

    vwin = get_browser_for_role(role, path);
    if (vwin != NULL) {
	gtk_window_present(GTK_WINDOW(vwin->main));
	return;
    }

    if (role == PKG_REGISTRY && n_user_handled_packages() == 0) {
	infobox(_("The gui package registry is empty"));
	return;
    }

    if (role == FUNC_FILES || role == NATIVE_DB) {
	title = files_title(role);
    } else if (role == PS_FILES) {
	title = g_strdup(_("gretl: example scripts"));
    } else if (role == TEXTBOOK_DATA) {
	title = g_strdup(_("gretl: data files"));
    } else if (role == REMOTE_DB) {
	title = g_strdup(_("gretl: databases on server"));
    } else if (role == DBNOMICS_TOP) {
	title = g_strdup(_("gretl: DB.NOMICS providers"));
    } else if (role == REMOTE_FUNC_FILES) {
	title = g_strdup(_("gretl: function packages on server"));
    } else if (role == REMOTE_DATA_PKGS) {
	title = g_strdup(_("gretl: data packages on server"));
    } else if (role == ADDONS_FILES) {
	title = g_strdup(_("gretl: addons"));
    } else if (role == PKG_REGISTRY) {
	title = g_strdup(_("gretl: packages on menus"));
    } else if (role == DBNOMICS_DB) {
	mypath = g_strdup(path);
	title = g_strdup_printf("gretl: %s datasets", path);
    } else if (role == DBNOMICS_SERIES) {
	mypath = g_strdup(path);
	title = g_strdup_printf("gretl: %s", path);
    }

    vwin = gretl_browser_new(role, title);
    g_free(title);

    if (role == REMOTE_DB) {
	gtk_window_set_default_size(GTK_WINDOW(vwin->main), 640, 480);
    }

    /* vertical box to hold file-listing widget and other elements */
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    make_files_toolbar(vwin);

    if (notebook_needed(role)) {
	/* we'll need more than one tab */
	filebox = files_notebook(vwin, role);
    } else {
	/* no tabs needed */
	filebox = files_vbox(vwin);
    }

    if (filebox == NULL) {
	gtk_widget_destroy(vwin->main);
	return;
    }

    gtk_box_pack_start(GTK_BOX(vwin->vbox), filebox, TRUE, TRUE, 0);

    if (role == TEXTBOOK_DATA) {
	file_collection *collection;
	GList *L = collections_for_role(role);

	build_datafiles_popup(vwin);
	while (L) {
	    collection = L->data;
	    g_signal_connect(G_OBJECT(collection->listbox), "button-press-event",
			     G_CALLBACK(popup_menu_handler),
			     vwin->popup);
	    L = L->next;
	}
    } else if (role == FUNC_FILES || role == REMOTE_FUNC_FILES ||
	       role == ADDONS_FILES || role == PKG_REGISTRY) {
	g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
			 G_CALLBACK(funcfiles_popup_handler),
			 vwin);
	if (role == FUNC_FILES) {
	    connect_menu_adjust_signal(vwin);
	}
    } else if (role == NATIVE_DB || role == REMOTE_DB) {
	build_db_popup(vwin);
	g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
			 G_CALLBACK(popup_menu_handler),
			 vwin->popup);
    } else if (role == REMOTE_DATA_PKGS)  {
	build_data_pkg_popup(vwin);
	g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
			 G_CALLBACK(popup_menu_handler),
			 vwin->popup);
    }

    if (REMOTE_ACTION(role)) {
	GtkWidget *hbox;

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
	if (DBNOMICS_ACTION(role)) {
	    vwin->status = gtk_label_new("");
	} else {
	    vwin->status = gtk_label_new(_("Network status: OK"));
	}
	gtk_label_set_justify(GTK_LABEL(vwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->status, FALSE, FALSE, 0);
    }

    /* put stuff into list box(es) */
    if (notebook_needed(role)) {
	err = populate_notebook_filelists(vwin, filebox, role);
    } else if (role == FUNC_FILES) {
	err = populate_filelist(vwin, NULL);
    } else if (role == NATIVE_DB) {
	gint w, h, ndb = 0;

	err = populate_filelist(vwin, &ndb);
	if (!err && ndb > 12) {
	    gtk_widget_get_size_request(filebox, &w, &h);
	    h += 100;
	    gtk_widget_set_size_request(filebox, w, h);
	}
    } else {
	err = populate_filelist(vwin, mypath);
    }

    if (err) {
	gtk_widget_destroy(vwin->main);
    } else {
	gtk_widget_show_all(vwin->main);
	gtk_widget_grab_focus(vwin->listbox);
	if (role == NATIVE_DB || role == FUNC_FILES) {
	    set_up_viewer_drag_target(vwin);
	}
	if (want_keyvals(role)) {
	    set_up_browser_keystrokes(vwin);
	}
    }

    if (err) {
	return;
    }

    if (role != TEXTBOOK_DATA && role != PS_FILES) {
	listbox_select_first(vwin);
    }
}

static int display_files_code (const gchar *s)
{
    if (!strcmp(s, "DisplayDataFiles"))
	return TEXTBOOK_DATA;
    if (!strcmp(s, "DisplayScripts"))
	return PS_FILES;
    if (!strcmp(s, "NativeDB"))
	return NATIVE_DB;
    if (!strcmp(s, "RemoteDB"))
	return REMOTE_DB;
    if (!strcmp(s, "LocalGfn"))
	return FUNC_FILES;
    if (!strcmp(s, "RemoteGfn"))
	return REMOTE_FUNC_FILES;
    if (!strcmp(s, "Addons"))
	return ADDONS_FILES;
    if (!strcmp(s, "DBNbrowse"))
	return DBNOMICS_TOP;

    return 0;
}

/* make a browser window to display a set of files: textbook
   data files, example (formerly: practice) scripts, databases...
*/

void show_files (GtkAction *action, gpointer p)
{
    int code = display_files_code(gtk_action_get_name(action));

    display_files(code, NULL);
}

void show_native_dbs (void)
{
    display_files(NATIVE_DB, NULL);
}

/* functions pertaining to gfn function packages */

static int populate_gfn_registry_list (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter;
    char *name;
    char *path;
    char *label;
    int modelwin;
    int i, n;

    store = GTK_LIST_STORE(gtk_tree_view_get_model
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    n = n_registered_packages();

    for (i=0; i<n; i++) {
	get_registered_pkg_info(i, &name, &path, &label, &modelwin);
	if (name != NULL && path != NULL) {
	    gchar *fullpath, *upath, *s = path;

	    if (!strncmp(s, "/menubar/", 9)) {
		s += 9;
	    }
	    upath = user_friendly_menu_path(s, modelwin);
	    if (upath != NULL) {
		fullpath = g_strdup_printf("%s/%s", upath, label);
	    } else {
		fullpath = g_strdup_printf("%s/%s", s, label);
	    }
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter,
			       0, name,
			       1, modelwin ? "Model" : "Main",
			       2, fullpath,
			       -1);
	    g_free(upath);
	    g_free(fullpath);
	}
    }

    return 0;
}

static int get_func_info (const char *path, char **pdesc,
			  char **pver, char **pdate,
			  char **pauthor, int *pdfdoc)
{
    int err;

    err = get_function_file_header(path, pdesc, pver, pdate,
				   pauthor, pdfdoc);
    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static int is_functions_dir (const char *path)
{
    int n = strlen(path) - 9;

    return n > 0 && !strcmp(path + n, "functions");
}

/* For a function package that lives in its own directory,
   see if it has an "examples" subdir.
*/

static int have_examples (const char *dirname)
{
    struct stat sbuf;
    gchar *test;
    int ret = 0;

    test = g_strdup_printf("%s%cexamples", dirname, SLASH);

    if (stat(test, &sbuf) == 0 &&
	(sbuf.st_mode & S_IFDIR)) {
	ret = 1;
    }

    g_free(test);

    return ret;
}

char *maybe_ellipsize_string (char *s, int maxlen)
{
    size_t n = g_utf8_strlen(s, -1);

    if (n > maxlen) {
	gretl_utf8_truncate(s, maxlen - 3);
	strcat(s, "...");
    }

    return s;
}

static void browser_insert_gfn_info (const char *pkgname,
				     const char *version,
				     const char *date,
				     const char *author,
				     const char *summary,
				     const char *dirname,
				     int uses_subdir,
				     int pdfdoc,
				     GtkListStore *store,
				     GtkTreeIter *iter)
{
    gint flags = 0;

    if (uses_subdir) {
	if (have_examples(dirname)) {
	    flags |= PKG_ATTR_RES;
	}
	if (pdfdoc) {
	    flags |= PKG_ATTR_DOC;
	}
    }

    gtk_list_store_set(store, iter,
		       0, pkgname,
		       1, version,
		       2, date,
		       3, author,
		       4, summary,
		       5, dirname,
		       6, flags,
		       -1);
}

static void check_loaded_gfn (const char *pkgname,
			      const char *fullname)
{
    const char *path;

    path = get_function_package_path_by_name(pkgname);

    if (path != NULL && strcmp(path, fullname)) {
	/* Watch out, there's a package of the same name already
	   loaded from a different location.
	*/
	fnpkg *pkg = get_function_package_by_filename(path, NULL);
	GtkWidget *editor = function_package_get_editor(pkg);

	if (editor != NULL) {
	    /* warn, don't unload the package being edited */
	    gchar *msg;

	    msg = g_strdup_printf("%s: package in new browser window may\n"
				  "conflict with package editor", pkgname);
	    gtk_window_present(GTK_WINDOW(editor));
	    msgbox(msg, GTK_MESSAGE_WARNING, editor);
	    g_free(msg);
	} else {
	    function_package_unload_full_by_filename(path);
	}
    }
}

/* Run various checks on @fullname: if all is OK, return 1; if
   something is amiss, return 0.
*/

static int ok_gfn_path (const char *fullname,
			const char *shortname,
			const char *dirname,
			GtkListStore *store,
			GtkTreeIter *iter,
			int subdir)
{
    char *descrip = NULL;
    char *version = NULL;
    char *date = NULL;
    char *author = NULL;
    int pdfdoc = 0;
    int err, ok = 0;

    /* Note that even if this is a dry run with @store = NULL,
       it may be worth performing the next action as a sanity
       check on the purported gfn.
    */
    err = get_func_info(fullname, &descrip, &version, &date,
			&author, &pdfdoc);

#if GFN_DEBUG > 1
    fprintf(stderr, "%s: %s: err=%d\n", dirname, shortname, err);
#endif

    if (!err) {
	if (store != NULL && iter != NULL) {
	    /* actually enter the file into the browser */
	    gchar *pkgname = g_strndup(shortname, strlen(shortname) - 4);

	    if (gfn_is_loaded(shortname)) {
		check_loaded_gfn(pkgname, fullname);
	    }
	    gtk_list_store_append(store, iter);
	    browser_insert_gfn_info(pkgname,
				    version,
				    date,
				    author,
				    descrip,
				    dirname,
				    subdir,
				    pdfdoc,
				    store,
				    iter);
	    g_free(pkgname);
	}
	ok = 1;
    }

    free(descrip);
    free(version);
    free(date);
    free(author);

    return ok;
}

/* Read (or simply just count) the .gfn files in a given directory.
   The signal to count rather than read is that the @store argument
   is NULL.
*/

static void
read_fn_files_in_dir (GDir *dir, const char *path,
		      GtkListStore *store,
		      GtkTreeIter *iter,
		      int *nfn)
{
    const gchar *basename;
    char fullname[MAXLEN];

    /* Look first for a gfn file in its own subdir, as
       in functions/foo/foo.gfn. That way if a package
       has been updated to zipfile status and there's
       also an older "plain gfn" version lying around
       we should get the newer one.
    */

    if (is_functions_dir(path)) {
	while ((basename = g_dir_read_name(dir)) != NULL) {
	    if (!strcmp(basename, ".") ||
		!strcmp(basename, "..")) {
		continue;
	    }
	    gretl_build_path(fullname, path, basename, NULL);
	    if (gretl_isdir(fullname)) {
		/* construct functions/foo/foo.gfn */
		gchar *realbase, *realpath;

		strcat(fullname, SLASHSTR);
		strcat(fullname, basename);
		strcat(fullname, ".gfn");
		if (gretl_file_exists(fullname)) {
		    realbase = g_strdup_printf("%s.gfn", basename);
		    realpath = g_strdup_printf("%s%c%s", path, SLASH, basename);
		    *nfn += ok_gfn_path(fullname, realbase, realpath,
					store, iter, 1);
		    g_free(realbase);
		    g_free(realpath);
		} else {
		    gretl_error_clear();
		}
	    }
	}
	g_dir_rewind(dir);
    }

    /* then look for "plain gfn" files */

    while ((basename = g_dir_read_name(dir)) != NULL) {
	if (!strcmp(basename, ".") ||
	    !strcmp(basename, "..")) {
	    continue;
	}
	if (has_suffix(basename, ".gfn")) {
	    gretl_build_path(fullname, path, basename, NULL);
	    *nfn += ok_gfn_path(fullname, basename, path,
				store, iter, 0);
	}
    }
}

#if GFN_DEBUG

static void show_dirs_list (char **S, int n, const char *msg)
{
    int i;

    fprintf(stderr, "*** dirs list: %s\n", msg);
    for (i=0; i<n; i++) {
	fprintf(stderr, " %d: '%s'\n", i, S[i]);
    }
}

#endif

/* When two gfns with the same name, @s, have been found, determine
   if one or other should be deleted, or perhaps just ignored.
   Our handles to the two packages are given by the GtkTreeIters
   @ia and @ib.
*/

static void maybe_delete_gfn_duplicate (GtkTreeModel *model,
					GtkTreeIter *ia,
					GtkTreeIter *ib,
					const char *s)

{
    gchar *va, *vb;
    gchar *da, *db;
    double vdiff;
    int del = 0;

    /* get version and date info for the two */
    gtk_tree_model_get(model, ia, 1, &va, 2, &da, -1);
    gtk_tree_model_get(model, ib, 1, &vb, 2, &db, -1);

    vdiff = dot_atof(va) - dot_atof(vb);

    if (vdiff > 0) {
	/* the second is older by version, delete it */
	del = 2;
    } else if (vdiff < 0) {
	/* the first is older by version, delete it */
	del = 1;
    } else {
	guint32 ed1 = get_epoch_day(da);
	guint32 ed2 = get_epoch_day(db);
	int ddiff = (int) ed1 - (int) ed2;

	if (ddiff > 0) {
	    /* the second is older by date, delete it */
	    del = 2;
	} else if (ddiff < 0) {
	    /* the first is older by date, delete it */
	    del = 1;
	} else {
	    /* same version and date: let's ignore, but not
	       actually delete, the second (per-user) version
	    */
	    del = -2;
	}
    }

    if (del != 0) {
	GtkListStore *store = GTK_LIST_STORE(model);
	GtkTreeIter *iter = del == 1 ? ia : ib;
	gchar *dirname = NULL;
	gchar *delpath = NULL;

	gtk_tree_model_get(model, iter, GFN_DIRNAME_COL, &dirname, -1);
	if (del > 0) {
	    delpath = g_strdup_printf("%s%c%s.gfn", dirname, SLASH, s);
	    gretl_remove(delpath);
	}
	gtk_list_store_remove(store, iter);
	g_free(dirname);
	g_free(delpath);
    }

    g_free(va);
    g_free(vb);
    g_free(da);
    g_free(db);
}

static void prune_gfn_duplicates (windata_t *vwin)
{
    GtkTreeModel *model;
    GtkTreeIter ia, ib;
    gchar *sa, *sb;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
    gtk_tree_model_get_iter_first(model, &ia);
    ib = ia;

    while (gtk_tree_model_iter_next(model, &ib)) {
	gtk_tree_model_get(model, &ia, 0, &sa, -1);
	gtk_tree_model_get(model, &ib, 0, &sb, -1);
	if (!strcmp(sa, sb)) {
	    maybe_delete_gfn_duplicate(model, &ia, &ib, sa);
	}
	g_free(sa);
	g_free(sb);
	ia = ib;
    }
}

/* Populate browser displaying gfn files installed on local machine:
   this is always called with a "clean slate": either we're showing a
   new browser window, or we're recreating the default listing after
   the user has pointed the browser at another directory (in which
   case the prior listing has been cleared out by the time we get
   here).
*/

static gint populate_gfn_list (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter;
    char **dnames = NULL;
    int i, n_dirs = 0;
    int nfn = 0;
    int err = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* compile an array of names of directories to search */
    dnames = get_plausible_search_dirs(FUNCS_SEARCH, &n_dirs);

#if GFN_DEBUG
    show_dirs_list(dnames, n_dirs, "gfn search");
#endif

    for (i=0; i<n_dirs; i++) {
	GDir *dir = gretl_opendir(dnames[i]);

	if (dir != NULL) {
	    read_fn_files_in_dir(dir, dnames[i], store, &iter, &nfn);
	    g_dir_close(dir);
	}
    }

    /* we're done with the directory names */
    strings_array_free(dnames, n_dirs);

    if (nfn == 0) {
	/* we didn't find any gfn files */
	warnbox(_("No gretl function packages were found on this computer.\n"
		  "Please try /File/Functions packages/On server"));
	err = 1;
    } else {
	int dups = 0;

	presort_treelist(vwin, &dups);
	if (dups) {
	    /* we get here if a given package is found in both
	       the "system" and the "per-user" location
	    */
	    prune_gfn_duplicates(vwin);
	}
    }

    return err;
}

static int gfn_paths_match (const char *p0, const char *p1,
			    const char *pkgname)
{
    int ret = 0;

    if (!strcmp(p0, p1)) {
	ret = 1;
    } else {
	/* allow for the possibility that @p1 has had a
	   package-specific subdir appended, relative to @p0
	*/
	size_t n = strlen(p0);

	if (strlen(p1) > n && !strncmp(p1, p0, n) &&
	    IS_SLASH(p1[n])) {
	    ret = !strcmp(p1 + n + 1, pkgname);
	}
    }

    return ret;
}

void set_gfn_add_button_state (const char *pkgname,
                               windata_t *vwin,
                               gboolean state)
{
    gchar *name = NULL;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &name);
    if (name != NULL && !strcmp(name, pkgname)) {
        GtkWidget *button;

        button = g_object_get_data(G_OBJECT(vwin->mbar), "add-button");
        if (button != NULL) {
            gtk_widget_set_sensitive(button, state);
        }
    }
    g_free(name);
}

static void update_gfn_browser (const char *pkgname,
				const char *version,
				const char *date,
				const char *author,
				const char *descrip,
				const char *fname,
				int uses_subdir,
				int pdfdoc,
				windata_t *vwin)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *p, *dirname;
    gchar *summary;
    gchar *oldname, *oldver, *olddir;
    int dirmatch = 0;
    int done = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
    if (!gtk_tree_model_get_iter_first(model, &iter)) {
	return;
    }

    dirname = g_strdup(fname);
    p = strrslash(dirname);
    if (p != NULL) {
	*p = '\0';
    }

    /* in case of truncation */
    summary = g_strdup(descrip);

    while (1) {
	gtk_tree_model_get(model, &iter, 0, &oldname, 1, &oldver,
			   GFN_DIRNAME_COL, &olddir, -1);
	if (!strcmp(oldname, pkgname) && !strcmp(oldver, version)) {
	    /* Found a match for package name and version: so update
	       the browser entry and record that we're done.
	    */
	    fprintf(stderr, "gfn update: updating %s %s\n", pkgname, version);
	    browser_insert_gfn_info(pkgname, version, date, author, summary,
				    dirname, uses_subdir, pdfdoc,
				    GTK_LIST_STORE(model), &iter);
	    done = 1;
	} else if (!dirmatch) {
	    dirmatch = gfn_paths_match(olddir, dirname, pkgname);
	}
	g_free(oldname); g_free(oldver); g_free(olddir);
	if (done || !gtk_tree_model_iter_next(model, &iter)) {
	    break;
	}
    }

    if (!done && dirmatch) {
	/* We didn't find an entry that matched by pkgname and
	   version, but we did determine that the browser was
	   pointing at a directory in which the package in
	   question would be found, if it were re-read. So it
	   seems we should append the package (and re-sort the
	   package list).
	*/
	fprintf(stderr, "gfn update: appending %s %s\n", pkgname, version);
	gtk_list_store_append(GTK_LIST_STORE(model), &iter);
	browser_insert_gfn_info(pkgname, version, date, author, summary,
				dirname, uses_subdir, pdfdoc,
				GTK_LIST_STORE(model), &iter);
	presort_treelist(vwin, NULL);
    }

    g_free(dirname);
    g_free(summary);
}

/* Update function package status, if needed, either after
   a call to save a function package, or after deleting a
   package by CLI means; the latter case is flagged by
   NULL values for @version and @descrip.
*/

void maybe_update_gfn_browser (const char *pkgname,
			       const char *version,
			       const char *date,
			       const char *author,
			       const char *descrip,
			       const char *fname,
			       int uses_subdir,
			       int pdfdoc)
{
    windata_t *vwin = get_browser_for_role(FUNC_FILES, NULL);
    int del = (version == NULL && descrip == NULL);

    if (vwin != NULL && vwin->listbox != NULL) {
	if (del) {
	    browser_delete_row_by_content(vwin, 0, pkgname,
					  GFN_DIRNAME_COL, fname);
	} else {
	    update_gfn_browser(pkgname, version, date, author, descrip,
			       fname, uses_subdir, pdfdoc, vwin);
	}
    }
}

void maybe_update_pkg_registry_window (const char *pkgname,
				       int code)
{
    windata_t *vwin = get_browser_for_role(PKG_REGISTRY, NULL);

    if (vwin != NULL && vwin->listbox != NULL) {
	if (code == MENU_ADD_FN_PKG) {
	    populate_gfn_registry_list(vwin);
	} else if (code == MENU_REMOVE_FN_PKG) {
	    browser_delete_row_by_content(vwin, 0, pkgname,
					  0, NULL);
	}
    }
}

gint populate_filelist (windata_t *vwin, gpointer p)
{
    if (vwin->role == NATIVE_DB) {
	return populate_dbfilelist(vwin, p);
    } else if (vwin->role == REMOTE_DB) {
	return populate_remote_db_list(vwin);
    } else if (vwin->role == DBNOMICS_TOP) {
	return populate_dbnomics_provider_list(vwin);
    } else if (vwin->role == DBNOMICS_DB) {
	return populate_dbnomics_dataset_list(vwin, p);
    } else if (vwin->role == DBNOMICS_SERIES) {
	return populate_dbnomics_series_list(vwin, p);
    } else if (vwin->role == REMOTE_FUNC_FILES) {
	return populate_remote_func_list(vwin, 0);
    } else if (vwin->role == REMOTE_DATA_PKGS) {
	return populate_remote_data_pkg_list(vwin);
    } else if (vwin->role == FUNC_FILES) {
	return populate_gfn_list(vwin);
    } else if (vwin->role == ADDONS_FILES) {
	return populate_addons_list(vwin);
    } else if (vwin->role == PKG_REGISTRY) {
	return populate_gfn_registry_list(vwin);
    } else {
	return read_file_descriptions(vwin, p);
    }
}

static GtkWidget *files_vbox (windata_t *vwin)
{
    const char *data_titles[] = {
	N_("File"),
	N_("Summary"),
	N_("Type")
    };
    const char *remote_data_titles[] = {
	N_("File"),
	N_("Source"),
	N_("Date")
    };
    const char *ps_titles[] = {
	N_("Script"),
	N_("Topic"),
	N_("Data")
    };
    const char *db_titles[] = {
	N_("Database"),
	N_("Source")
    };
    const char *remote_db_titles[] = {
	N_("Database"),
	N_("Source"),
	N_("Local status")
    };
    const char *dbnomics_top_titles[] = {
	N_("Code"),
	N_("Name")
    };
    const char *dbnomics_db_titles[] = {
	N_("Code"),
	N_("Content")
    };
    const char *dbnomics_series_titles[] = {
	N_("Code"),
	N_("Description")
    };
    const char *func_titles[] = {
	N_("Package"),
	N_("Version"),
	N_("Date"),
	N_("Author"),
	N_("Summary")
    };
    const char *remote_func_titles[] = {
	N_("Package"),
	N_("Version"),
	N_("Date"),
	N_("Author"),
	N_("Summary"),
	N_("Local status")
    };
    const char *addons_titles[] = {
	N_("Addon"),
	N_("Summary"),
	N_("Date")
    };
    const char *registry_titles[] = {
	N_("Package"),
	N_("Window"),
	N_("Menu")
    };

    GType types_2[] = {
	G_TYPE_STRING,
	G_TYPE_STRING
    };
    GType types_3[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING
    };
    GType addons_types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING    /* hidden string: full path */
    };
    GType func_types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,   /* hidden string: directory */
	G_TYPE_INT       /* hidden flags: has examples dir? doc? */
    };
    GType remote_func_types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_BOOLEAN,  /* hidden boolean: zipfile? */
	G_TYPE_STRING    /* hidden string: dependencies */
    };
    const char **titles = data_titles;
    GType *types = types_2;
    int full_width = 580, file_height = 300;
    int hidden_cols = 0;
    int use_tree = 0;
    GtkWidget *vbox;
    int cols = 2;

    switch (vwin->role) {
    case TEXTBOOK_DATA:
	titles = data_titles;
	cols = 3;
	full_width = 600;
	break;
    case NATIVE_DB:
	titles = db_titles;
	cols = 3;
	hidden_cols = 1;
	break;
    case REMOTE_DB:
	titles = remote_db_titles;
	cols = 3;
	use_tree = 1;
	break;
    case DBNOMICS_TOP:
	titles = dbnomics_top_titles;
	full_width = 650;
	break;
    case DBNOMICS_DB:
	titles = dbnomics_db_titles;
	full_width = 650;
	break;
    case DBNOMICS_SERIES:
	titles = dbnomics_series_titles;
	full_width = 600;
	break;
    case REMOTE_DATA_PKGS:
	titles = remote_data_titles;
	cols = 3;
	full_width = 600;
	break;
    case PS_FILES:
	titles = ps_titles;
	cols = 3;
	full_width = 560;
	break;
    case FUNC_FILES:
	titles = func_titles;
	types = func_types;
	cols = G_N_ELEMENTS(func_types);
	hidden_cols = 2;
	full_width = 760;
	file_height = 320;
	break;
    case REMOTE_FUNC_FILES:
	titles = remote_func_titles;
	types = remote_func_types;
	cols = G_N_ELEMENTS(remote_func_types);
	hidden_cols = 2;
	full_width = 760;
	file_height = 340;
	break;
    case ADDONS_FILES:
	titles = addons_titles;
	types = addons_types;
	cols = G_N_ELEMENTS(addons_types);
	hidden_cols = 1;
	full_width = 520;
	break;
    case PKG_REGISTRY:
	titles = registry_titles;
	cols = 3;
	full_width = 600;
	break;
    default:
	break;
    }

    if (cols == 3) {
	types = types_3;
    }

    full_width *= gui_scale;
    file_height *= gui_scale;

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_widget_set_size_request(vbox, full_width, file_height);
    /* note: the following packs and attaches vwin->listbox */
    vwin_add_list_box(vwin, GTK_BOX(vbox), cols, hidden_cols,
		      types, titles, use_tree);
    gtk_widget_show(vbox);

    return vbox;
}

static void switch_files_page (GtkNotebook *notebook,
			       GtkWidget *page,
			       guint pgnum,
			       windata_t *vwin)
{
    GtkWidget *tab = gtk_notebook_get_nth_page(notebook, pgnum);

    vwin->listbox = g_object_get_data(G_OBJECT(tab), "listbox");
}

/* below: construct a set of notebook pages for either data file
   collections (Ramanathan, Wooldridge, etc.) or example scripts.
   The function creates the pages but does not yet fill them out.
*/

static GtkWidget *files_notebook (windata_t *vwin, int role)
{
    GList *L = NULL;
    file_collection *collection;
    GtkWidget *notebook;
    GtkWidget *page;
    GtkWidget *label;
    int err = 0;

    /* assemble the info we'll need */
    err = build_file_collections(role);
    if (err) {
	fprintf(stderr, "files_notebook: build_files_collections failed\n");
	return NULL;
    }

    L = collections_for_role(role);
    if (L == NULL) {
	return NULL;
    }

    notebook = gtk_notebook_new();
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook), TRUE);

    while (L) {
	collection = L->data;
	page = files_vbox(vwin);
	label = gtk_label_new(collection->title);
	gtk_widget_show(label);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
	collection->listbox = vwin->listbox;
	g_object_set_data(G_OBJECT(collection->listbox), "collection",
			  collection);
	g_object_set_data(G_OBJECT(page), "listbox", collection->listbox);
	g_signal_connect(G_OBJECT(collection->listbox), "key-press-event",
			 G_CALLBACK(enter_opens_file), vwin);
	L = L->next;
    }

    g_signal_connect(G_OBJECT(notebook), "switch-page",
		     G_CALLBACK(switch_files_page),
		     vwin);
    if (gtk_notebook_get_n_pages(GTK_NOTEBOOK(notebook)) > 5) {
	gtk_notebook_popup_enable(GTK_NOTEBOOK(notebook));
    }

    gtk_widget_show(notebook);

    return notebook;
}

/* below: fill out a set of notebook pages (for data files
   or script files), entering the details into the page,
   then select the page to display
*/

static int populate_notebook_filelists (windata_t *vwin,
					GtkWidget *notebook,
					int role)
{
    file_collection *collection;
    file_collection *selected = NULL;
    const char *page = NULL;
    GList *L = NULL;
    int found = 0;
    int pgnum = 0;

    L = collections_for_role(role);
    if (L == NULL) {
	return 1;
    }

    if (role == TEXTBOOK_DATA) {
	page = get_datapage();
    } else if (role == PS_FILES) {
	page = get_scriptpage();
    }

    selected = L->data;

    while (L) {
	collection = L->data;
	vwin->listbox = collection->listbox;
	populate_filelist(vwin, collection);
	if (page != NULL && !strcmp(collection->title, page)) {
	    selected = collection;
	    pgnum = found;
	}
	found++;
	L = L->next;
    }

    if (found == 0) {
	/* didn't find anything to show */
	fprintf(stderr, "populate_notebook_filelists: found = 0!\n");
	return 1;
    }

    vwin->listbox = selected->listbox;
    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), pgnum);
    gtk_widget_grab_focus(vwin->listbox);

    return 0;
}
