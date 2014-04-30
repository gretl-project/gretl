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

#include "gretl.h"
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

#include "gretl_xml.h"
#include "gretl_func.h"

#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>

static GtkWidget *files_vbox (windata_t *vwin);
static GtkWidget *files_notebook (windata_t *vwin, int role);
static int populate_notebook_filelists (windata_t *vwin, 
					GtkWidget *notebook,
					int role);

typedef struct _file_collection file_collection;

struct _file_collection {
    char *path;
    char *descfile;
    char *title;
    int which;
    GtkWidget *listbox;
};

enum {
    STACK_PUSH,
    STACK_POP_DATA,
    STACK_POP_PS,
    STACK_RESET_DATA,
    STACK_RESET_PS,
    STACK_SORT_DATA,
    STACK_SORT_PS,
    STACK_DESTROY
};

enum {
    COLL_NONE,
    COLL_DATA,
    COLL_PS
};

struct fpkg_response {
    int col1_width;
    int try_server;
};

#define REMOTE_ACTION(c) (c == REMOTE_DB || \
                          c == REMOTE_FUNC_FILES || \
                          c == REMOTE_DATA_PKGS || \
                          c == REMOTE_ADDONS)

static void
read_fn_files_in_dir (DIR *dir, const char *path, 
		      GtkListStore *store, GtkTreeIter *iter,
		      int *nfn, int *maxlen);

static void 
fpkg_response_init (struct fpkg_response *f, gpointer data)
{
    f->col1_width = 0;

    if (data == NULL || data == mdata) {
	/* called from main menu */
	f->try_server = 0;
    } else {
	f->try_server = -1;
    }
}

static char *full_path (char *s1, const char *s2)
{
    static char fpath[FILENAME_MAX];
    int n = strlen(s1);

    if (s1[n-1] == '.') {
	s1[n-1] = '\0';
	n--;
    }
    
    if (s1[n-1] == SLASH) {
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
		coll->which = COLL_DATA;
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
	    coll->which = COLL_PS;
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
    fp = gretl_fopen(test, "r");

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

    coll->which = COLL_NONE;
    coll->title = NULL;
    coll->path = gretl_strdup(path);
    coll->descfile = gretl_strdup(descfile);

    if (coll->path == NULL || coll->descfile == NULL) {
	*err = E_ALLOC;
    } else {
	int os = is_oldstyle_collection(coll, err);

	if (!*err && !os) {
	    if (strstr(coll->descfile, "ps_")) {
		coll->which = COLL_PS;
	    } else {
		coll->which = COLL_DATA;
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
    const file_collection *ca = *(const file_collection **) a;
    const file_collection *cb = *(const file_collection **) b;

    if (!strcmp(ca->title, "Gretl")) {
	return -1;
    } else if (!strcmp(cb->title, "Gretl")) {
	return 1;
    } else {
	return strcmp(ca->title, cb->title);
    }
}

static void collection_stack_sort (file_collection **colls, int n)
{
    if (n >= 2) {
	qsort(colls, n, sizeof *colls, compare_colls);
    }
}

static file_collection *collection_stack (file_collection *coll, int op)
{
    static file_collection **datacoll;
    static file_collection **pscoll;
    static int n_data;
    static int n_data_popped;
    static int n_ps;
    static int n_ps_popped;
    file_collection *ret = NULL;
    int j;

    if (op == STACK_PUSH && coll != NULL) {
	if (coll->which == COLL_DATA) {
	    datacoll = realloc(datacoll, (n_data + 1) * sizeof *datacoll);
	    if (datacoll != NULL) {
		datacoll[n_data++] = coll;
		ret = coll;
	    }
	} else if (coll->which == COLL_PS) {
	    pscoll = realloc(pscoll, (n_ps + 1) * sizeof *pscoll);
	    if (pscoll != NULL) {
		pscoll[n_ps++] = coll;
		ret = coll;
	    }
	}
    } else if (op == STACK_POP_DATA && n_data_popped < n_data) {
	ret = datacoll[n_data_popped++];
    } else if (op == STACK_POP_PS && n_ps_popped < n_ps) {
	ret = pscoll[n_ps_popped++];
    } else if (op == STACK_RESET_DATA) {
	n_data_popped = 0;
    } else if (op == STACK_RESET_PS) {
	n_ps_popped = 0;
    } else if (op == STACK_SORT_DATA) {
	collection_stack_sort(datacoll, n_data);
    } else if (op == STACK_SORT_PS) {
	collection_stack_sort(pscoll, n_ps);
    } else if (op == STACK_DESTROY) {
        for (j=0; j<n_data; j++) {
	    free_file_collection(datacoll[j]);
	}
        free(datacoll);
        datacoll = NULL;
        n_data = 0;
        n_data_popped = 0;

        for (j=0; j<n_ps; j++) {
	    free_file_collection(pscoll[j]);
	}
        free(pscoll);
        pscoll = NULL;
        n_ps = 0;
        n_ps_popped = 0;
    } 

    return ret;
}

static int push_collection (file_collection *collection)
{
    return (collection_stack(collection, STACK_PUSH) == NULL);
}

static file_collection *pop_file_collection (int role)
{
    if (role == TEXTBOOK_DATA) {
	return collection_stack(NULL, STACK_POP_DATA);
    } else {
	return collection_stack(NULL, STACK_POP_PS);
    }
}

void destroy_file_collections (void)
{
    collection_stack(NULL, STACK_DESTROY);
    maybe_add_packages_to_model_menus(NULL);
}

static void reset_files_stack (int role)
{
    if (role == TEXTBOOK_DATA) {
	collection_stack(NULL, STACK_RESET_DATA);
    } else {
	collection_stack(NULL, STACK_RESET_PS);
    }
}

static void sort_files_stack (int role)
{
    if (role == TEXTBOOK_DATA) {    
	collection_stack(NULL, STACK_SORT_DATA);
    } else {
	collection_stack(NULL, STACK_SORT_PS);
    }
}

/* Returns the number of file collections found and pushed;
   writes non-zero to @err if something show-stopping
   occurs
*/

static int get_file_collections_from_dir (const char *path, DIR *dir,
					  int *err)
{
    file_collection *coll;
    struct dirent *dirent;
    int n = 0;

    while (!*err && (dirent = readdir(dir))) { 
	/* we're looking for a filename that ends with "descriptions" */
	if (strstr(dirent->d_name, "descriptions")) {
	    size_t len = strlen(dirent->d_name);

#if COLL_DEBUG
	    fprintf(stderr, "   %s: looking at '%s'\n", path, dirent->d_name);
#endif
	    if (!strcmp(dirent->d_name + len - 12, "descriptions")) {
		coll = file_collection_new(path, dirent->d_name, err);
		if (coll != NULL) {
		    *err = push_collection(coll);
		    if (!*err) {
			n++;
		    }
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
    char *path = NULL;
    DIR *topdir;      
    struct dirent *dirent;
    int n_coll = 0;

#if COLL_DEBUG
    fprintf(stderr, "*** seek_file_collections: basedir='%s', type=%d\n", 
	    basedir, stype);
#endif

    if (*err) {
	/* a fatal error occurred already, skip it */
	return 0;
    }

    if (stype == DATA_SEARCH) {
	path = gretl_strdup_printf("%sdata", basedir);
    } else if (stype == SCRIPT_SEARCH) {
	path = gretl_strdup_printf("%sscripts", basedir);
    } else {
	/* USER_SEARCH */
	path = gretl_strdup(basedir);
	trim_slash(path);
    } 

    topdir = gretl_opendir(path);
    if (topdir == NULL) {
	free(path);
	return 0;
    }

#if COLL_DEBUG
    fprintf(stderr, "*** seek_file_collections: path='%s'\n", path);
#endif

    while (!*err && (dirent = readdir(topdir))) {
	if (!dont_go_there(dirent->d_name)) {
	    char *subpath;
	    DIR *subdir;

#if COLL_DEBUG > 1
	    fprintf(stderr, " dname = '%s'\n", dirent->d_name);
#endif
	    if (strcmp(dirent->d_name, ".")) {
		subpath = full_path(path, dirent->d_name);
	    } else {
		subpath = path;
	    }
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
		closedir(subdir);
	    }
	}
    }

    closedir(topdir);
    free(path);

#if COLL_DEBUG
    fprintf(stderr, "*** found %d collections\n", n_coll);
#endif

    return n_coll;
}

#if COLL_DEBUG
static void print_collection (const file_collection *coll)
{
    printf("path = '%s'\n", coll->path);
    printf("descfile = '%s'\n", coll->descfile);
    if (coll->title != NULL && *coll->title != '\0') {
	printf("title = '%s'\n", coll->title);
    }
}

static void print_collections (int role)
{
    file_collection *coll;

    if (role == TEXTBOOK_DATA) {
	printf("\n*** Data collections:\n");
    } else {
	printf("\n*** Script collections:\n");
    }

    while ((coll = pop_file_collection(role))) {
	print_collection(coll);
    }

    reset_files_stack(role);
}
#endif

static int build_file_collections (void)
{
    static int built;
    static int err;

    if (!built && !err) {
	const char *wd;
	int n = 0;

	n += seek_file_collections(gretl_home(), DATA_SEARCH, &err);
	n += seek_file_collections(gretl_home(), SCRIPT_SEARCH, &err);
#ifdef OS_OSX
	n += seek_file_collections(gretl_app_support_dir(), DATA_SEARCH, &err);
	n += seek_file_collections(gretl_app_support_dir(), SCRIPT_SEARCH, &err);
#else
	n += seek_file_collections(gretl_dotdir(), DATA_SEARCH, &err);
	n += seek_file_collections(gretl_dotdir(), SCRIPT_SEARCH, &err);
#endif
	n += seek_file_collections(gretl_workdir(), USER_SEARCH, &err);
	if (!err) {
	    wd = maybe_get_default_workdir();
	    if (wd != NULL) {
		n += seek_file_collections(wd, USER_SEARCH, &err);
		n += seek_file_collections(wd, DATA_SEARCH, &err);
		n += seek_file_collections(wd, SCRIPT_SEARCH, &err);
	    }
	}
	if (!err && n == 0) {
	    err = E_DATA;
	}
	if (!err) {
	    sort_files_stack(TEXTBOOK_DATA);
	    sort_files_stack(PS_FILES);
	}
	built = 1;
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

	if (*line == '#') continue;

	if (win->role == TEXTBOOK_DATA) {
	    if (sscanf(line, " \"%23[^\"]\",\"%79[^\"]\"", 
		       fname, descrip) == 2) {
		err = validate_desc_strings(fname, descrip, NULL);
		if (!err) {
		    gtk_list_store_append(store, &iter);
		    gtk_list_store_set(store, &iter, 
				       0, strip_extension(fname), 
				       1, descrip, -1);
		}
	    }
	} else { 
	    /* script files */
	    if (sscanf(line, " \"%23[^\"]\",\"%79[^\"]\",\"%63[^\"]\"", 
		       fname, descrip, data) == 3) {
		err = validate_desc_strings(fname, descrip, data);
		if (!err) {
		    gtk_list_store_append(store, &iter);
		    gtk_list_store_set(store, &iter, 
				       0, strip_extension(fname), 
				       1, descrip, 
				       2, data, -1);
		}
	    }
	}
    }

    fclose(fp);

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

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &filename);
    collection = g_object_get_data(G_OBJECT(vwin->listbox), "collection");
    build_path(fullname, collection->path, filename, ".gdt");
    g_free(filename);

#if 0
    fprintf(stderr, "info: active=%d, fullname='%s'\n", vwin->active_var,
	    fullname);
    fprintf(stderr, "collection path='%s'\n", collection->path);
#endif

    descrip = gretl_get_gdt_description(fullname);

    if (descrip == NULL) {
	errbox(_("Failed to retrieve description of data"));
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
    gchar *filename;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &filename);
    collection = g_object_get_data(G_OBJECT(vwin->listbox), "collection");
    build_path(tryfile, collection->path, filename, ".gdt");
    g_free(filename);

    set_datapage(collection->title);

    verify_open_data(vwin, OPEN_DATA);
}

void browser_open_ps (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    file_collection *collection;
    gchar *filename;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &filename);
    collection = g_object_get_data(G_OBJECT(vwin->listbox), "collection");
    build_path(scriptfile, collection->path, filename, ".inp");
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

static int gui_load_user_functions (const char *fname)
{
    int err = load_function_package_by_filename(fname, NULL);

    if (err) {
	gui_errmsg(err);
    } 

    return err;
}

static int gui_delete_fn_pkg (const char *pkgname, const char *fname, 
			      windata_t *vwin)
{
    char *msg = g_strdup_printf(_("Function package %s"), fname);
    const char *opts[] = {
	N_("Delete package file"),
	N_("Unload member functions")
    };
    int active[] = {1, 1};
    int resp, err = 0;

    resp = checks_only_dialog("gretl", msg, opts, 2, 
			      active, 0, vwin->main);
    g_free(msg);

    if (resp < 0 || (active[0] == 0 && active[1] == 0)) {
	/* canceled, or equivalent */
	return 0;
    }

    /* remove entry from packages.xml, if present */
    unregister_function_package(pkgname);

    if (active[1]) {
	/* unload the package and its members from memory */
	function_package_unload_full_by_filename(fname);
    } else {
	/* unload the package (only) from memory */
	function_package_unload_by_filename(fname);
    }

     if (active[0]) {
	/* scratch the package file */
	err = gretl_remove(fname);
	if (err) {
	    file_write_errbox(fname);
	}
    }

    if (!err) {
	/* remove package from GUI listing */
	GtkTreeModel *mod;
	GtkTreeIter iter;
	int i = 0;

	mod = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
	gtk_tree_model_get_iter_first(mod, &iter);
	while (i < vwin->active_var) {
	    gtk_tree_model_iter_next(mod, &iter);
	    i++;
	}
	gtk_list_store_remove(GTK_LIST_STORE(mod), &iter);
    }

    return err;
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
	if (g_path_is_absolute(path) && !strstr(path, "dltmp.")) {
	    pprintf(prn, "File: %s\n", path);
	}
	err = print_function_package_info(path, prn);
    } else {
	err = print_function_package_code(path, prn);
    }
	
    if (err) {
	gretl_print_destroy(prn);
	gui_errmsg(err);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", pkgname);

	vwin = view_buffer(prn, 78, 350, title, role, NULL);
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

void set_funcs_dir_callback (windata_t *vwin, char *path)
{
    DIR *dir = gretl_opendir(path);

    if (dir != NULL) {
	const char *opts[] = {
	    N_("Add to functions shown"),
	    N_("Replace functions shown")
	};
	GtkListStore *store;
	GtkTreeIter iter;
	int maxlen = 0;
	int nfn = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(vwin->listbox), "nfn"));
	int nfn0 = nfn;	
	int resp;

	store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));

	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	read_fn_files_in_dir(dir, path, store, NULL, &nfn, &maxlen);

	if (nfn == nfn0) {
	    warnbox(_("No additional function files were found"));
	    closedir(dir);
	    return;
	} 

	resp = radio_dialog("gretl", "function packages", opts, 2, 
			    0, 0, vwin->main);

	if (resp < 0) {
	    /* canceled */
	    closedir(dir);
	    return;
	} else if (resp > 0) {
	    gtk_list_store_clear(store);
	    nfn = nfn0 = 0;
	} else {
	    nfn = nfn0;
	}

	rewinddir(dir);

	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	read_fn_files_in_dir(dir, path, store, &iter, &nfn, &maxlen);

	if (nfn > nfn0) {
	    GtkTreeSelection *sel;

	    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
	    gtk_tree_selection_selected_foreach(sel, fix_selected_row, vwin);
	    infobox(_("Found %d function file(s)"), nfn - nfn0);
	    g_object_set_data(G_OBJECT(vwin->listbox), "nfn", GINT_TO_POINTER(nfn));
	    if (resp > 0) {
		g_object_set_data(G_OBJECT(vwin->listbox), "keepdir", GINT_TO_POINTER(1));
	    }
	} else {
	    warnbox(_("No function files were found"));
	}
	    
	closedir(dir);
    }
}

static void browser_functions_handler (windata_t *vwin, int task)
{
    char path[FILENAME_MAX];
    gchar *pkgname = NULL;
    gchar *dir;
    int dircol = 0;
    int err = 0;

    if (vwin->role == FUNC_FILES) {
	dircol = 4;
    } else if (vwin->role != REMOTE_FUNC_FILES) {
	dummy_call();
	return;
    }

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &pkgname);

    if (dircol != 0) {
	tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			     dircol, &dir);
	build_path(path, dir, pkgname, ".gfn");
	g_free(dir);
    } else {
	strcpy(path, pkgname);
    }

#if 0
    fprintf(stderr, "browser_functions_handler: active=%d, pkgname='%s'\n"
	    "path='%s'\n", vwin->active_var, pkgname, path);
#endif

    if (task == LOAD_FN_PKG) {
	err = gui_load_user_functions(path);
    } else if (task == DELETE_FN_PKG) {
	err = gui_delete_fn_pkg(pkgname, path, vwin);
    } else if (task == VIEW_FN_PKG_INFO) {
	display_function_package_data(pkgname, path, VIEW_PKG_INFO);
    } else if (task == VIEW_FN_PKG_CODE) {
	display_function_package_data(pkgname, path, VIEW_PKG_CODE);
    } else if (task == EDIT_FN_PKG) {
	edit_function_package(path);
    } else if (task == MENU_ADD_FN_PKG) {
	gui_add_package_to_menu(path, TRUE, vwin->main);
    } else if (task == CALL_FN_PKG) {
	/* note: this is the double-click default */
	call_function_package(path, vwin, &err);
	if (err == FN_NO_DATA) {
	    display_function_package_data(pkgname, path, VIEW_PKG_INFO);
	    err = 0;
	}
    }

    g_free(pkgname);

    if (!err && task == CALL_FN_PKG) {
	maybe_update_func_files_window(task);
    }
}

static void show_addon_info (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *pkgname = NULL;
    gchar *descrip = NULL;
    gchar *status = NULL;
    int v = vwin->active_var;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), v, 0, &pkgname);
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), v, 4, &descrip);
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), v, 3, &status);
    
    if (pkgname == NULL || descrip == NULL) {
	gui_errmsg(E_DATA);
    } else {
	gchar *local = NULL;

	if (status != NULL && !strcmp(status, _("Not up to date"))) {
	    char *path = gretl_function_package_get_path(pkgname, PKG_SUBDIR);
	    gchar *ver = NULL, *date = NULL;
	    fnpkg *pkg = NULL;
	    int err = 0;

	    if (path != NULL) {
		pkg = get_function_package_by_filename(path, &err);
		err = function_package_get_properties(pkg, "version", &ver,
						      "date", &date,
						      NULL);
		if (!err) {
		    local = g_strdup_printf("Installed version is %s (%s)", 
					    ver, date);
		    g_free(ver);
		    g_free(date);
		}
		free(path);
	    }
	}

	if (local != NULL) {
	    infobox("%s:\n%s\n%s", pkgname, descrip, local);
	    g_free(local);
	} else {
	    infobox("%s:\n%s", pkgname, descrip);
	}
    }

    g_free(pkgname);
    g_free(descrip);
    g_free(status);
}

static void install_addon_callback (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *pkgname = NULL;
    int v = vwin->active_var;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), v, 0, &pkgname);
    
    if (pkgname == NULL) {
	gui_errmsg(E_DATA);
    } else {
	char *local_path = NULL;
	int err = download_addon(pkgname, &local_path);

	if (!err) {
	    list_store_set_string(GTK_TREE_VIEW(vwin->listbox), v, 3,
				  _("Up to date"));
	    /* if there was an old version of the addon loaded,
	       we need to unload it now so as to get correct 
	       information in response to a subsequent call to
	       "check for addons"
	    */
	    function_package_unload_full_by_filename(local_path);
	    free(local_path);
	} 

	g_free(pkgname);
    }
}

void browser_load_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, LOAD_FN_PKG);
} 

void browser_edit_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, EDIT_FN_PKG);
}

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

windata_t *get_local_viewer (int remote_role)
{
    windata_t *vwin = NULL;

    if (remote_role == REMOTE_DB) {
	vwin = get_browser_for_role(NATIVE_DB);
    } else if (remote_role == REMOTE_FUNC_FILES) {
	vwin = get_browser_for_role(FUNC_FILES);
    }

    return vwin;
}

static void new_package_callback (GtkWidget *w, gpointer p)
{
    functions_selection_wrapper();
}

static void close_files_viewer (GtkWidget *w, windata_t *vwin)
{
    gtk_widget_destroy(vwin->main);
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
    int dircol = 4;
    int ret = 0;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &pkgname);
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 dircol, &dirname);

    if (pkgname != NULL && dirname != NULL) {
	char path[FILENAME_MAX];

	build_path(path, dirname, pkgname, ".gfn");
	ret = package_is_available_for_menu(pkgname, path);
    }

    g_free(pkgname);
    g_free(dirname);

    return ret;
}

static void check_add_button_state (GtkTreeSelection *sel, windata_t *vwin)
{
    GtkWidget *button;

    button = g_object_get_data(G_OBJECT(vwin->mbar), "add-button");
    if (button != NULL) {
	gtk_widget_set_sensitive(button, get_menu_add_ok(vwin));
    }
}

static void connect_menu_add_signal (windata_t *vwin)
{
    GtkTreeSelection *sel;

    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));

    g_signal_connect(G_OBJECT(sel), "changed",
		     G_CALLBACK(check_add_button_state), 
		     vwin);
}

static void build_funcfiles_popup (windata_t *vwin)
{
    vwin->popup = gtk_menu_new();

    if (vwin->role == FUNC_FILES) {
	/* local function files: full menu */
	int menu_add_ok = 0;
	GtkWidget *add;

	add = g_object_get_data(G_OBJECT(vwin->mbar), "add-button");
	if (add != NULL && gtk_widget_is_sensitive(add)) {
	    menu_add_ok = 1;
	}

	add_popup_item(_("Edit"), vwin->popup, 
		       G_CALLBACK(browser_edit_func), 
		       vwin);
	add_popup_item(_("Info"), vwin->popup, 
		       G_CALLBACK(show_function_info), 
		       vwin);
	add_popup_item(_("View code"), vwin->popup, 
		       G_CALLBACK(show_function_code), 
		       vwin);
	add_popup_item(_("Execute"), vwin->popup, 
		       G_CALLBACK(browser_call_func), 
		       vwin);
	if (menu_add_ok) {
	    add_popup_item(_("Add to menu"), vwin->popup, 
			   G_CALLBACK(add_func_to_menu), 
			   vwin);
	}
	add_popup_item(_("Delete"), vwin->popup, 
		       G_CALLBACK(browser_del_func), 
		       vwin);
	add_popup_item(_("New"), vwin->popup, 
		       G_CALLBACK(new_package_callback), 
		       vwin);
    } else if (vwin->role == REMOTE_FUNC_FILES) {
	/* files on server: limited menu */
	add_popup_item(_("Info"), vwin->popup, 
		       G_CALLBACK(pkg_info_from_server), 
		       vwin);
	add_popup_item(_("Install"), vwin->popup, 
		       G_CALLBACK(install_file_from_server), 
		       vwin);
    } else if (vwin->role == REMOTE_ADDONS) {
	add_popup_item(_("Info"), vwin->popup, 
		       G_CALLBACK(show_addon_info), 
		       vwin);
	add_popup_item(_("Install"), vwin->popup, 
		       G_CALLBACK(install_addon_callback), 
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
    display_files(REMOTE_DB, p);
}

static void show_local_dbs (GtkWidget *w, gpointer p)
{
    display_files(NATIVE_DB, p);
}

static void show_server_funcs (GtkWidget *w, gpointer p)
{
    display_files(REMOTE_FUNC_FILES, p);
}

static void show_server_data_pkgs (GtkWidget *w, gpointer p)
{
    display_files(REMOTE_DATA_PKGS, p);
}

static void show_local_funcs (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, p);
}

static void alt_funcs_dir (GtkWidget *w, windata_t *vwin)
{
    file_selector_with_parent(SET_FDIR, FSEL_DATA_VWIN, vwin, 
                              vwin->main);
}

static void alt_db_dir (GtkWidget *w, windata_t *vwin)
{
    file_selector_with_parent(SET_DBDIR, FSEL_DATA_VWIN, vwin, 
                              vwin->main);
}

enum {
    BTN_EDIT = 1,
    BTN_INFO,
    BTN_CODE,
    BTN_INDX,
    BTN_INST,
    BTN_EXEC,
    BTN_ADD,
    BTN_DEL,
    BTN_WWW,
    BTN_HOME,
    BTN_NEW,
    BTN_FIND,
    BTN_OPEN,
    BTN_DIR,
    BTN_CLOSE
};

static GretlToolItem files_items[] = {
    { N_("Open"),           GTK_STOCK_OK, NULL, BTN_OPEN },
    { N_("Select directory"), GTK_STOCK_DIRECTORY, NULL, BTN_DIR },
    { N_("Edit"),           GTK_STOCK_EDIT,       G_CALLBACK(browser_edit_func), BTN_EDIT },
    { N_("Info"),           GTK_STOCK_INFO,       NULL,                          BTN_INFO },
    { N_("View code"),      GTK_STOCK_PROPERTIES, G_CALLBACK(show_function_code), BTN_CODE },
    { N_("List series"),    GTK_STOCK_INDEX,      NULL,                          BTN_INDX },
    { N_("Install"),        GTK_STOCK_SAVE,       NULL,                          BTN_INST },
    { N_("Execute"),        GTK_STOCK_EXECUTE,    G_CALLBACK(browser_call_func), BTN_EXEC },
    { N_("Add to menu"),    GTK_STOCK_ADD,        G_CALLBACK(add_func_to_menu),  BTN_ADD },
    { N_("Delete"),         GTK_STOCK_DELETE,     G_CALLBACK(browser_del_func),  BTN_DEL },
    { N_("Look on server"), GTK_STOCK_NETWORK,    NULL,                          BTN_WWW },
    { N_("Local machine"),  GTK_STOCK_HOME,       NULL,                          BTN_HOME },
    { N_("New"),            GTK_STOCK_NEW,        G_CALLBACK(new_package_callback), BTN_NEW },
    { N_("Windows"),        GRETL_STOCK_WINLIST,  GNULL, 0 },
    { N_("Close"),          GTK_STOCK_CLOSE,      G_CALLBACK(close_files_viewer), BTN_CLOSE }
};

static int n_files_items = G_N_ELEMENTS(files_items);

#define common_item(f) (f == 0 || f == BTN_FIND || f == BTN_CLOSE)

#define local_funcs_item(f) (f == BTN_EDIT || f == BTN_NEW || \
			     f == BTN_DEL || f == BTN_CODE)

static int files_item_get_callback (GretlToolItem *item, int role)
{
    if (common_item(item->flag)) {
	return 1;
    } else if (local_funcs_item(item->flag)) {
	return (role == FUNC_FILES);
    } else if (item->flag == BTN_INST) {
	if (role == REMOTE_ADDONS) {
	    item->func = G_CALLBACK(install_addon_callback);
	    return 1;
	} else {
	    item->func = G_CALLBACK(install_file_from_server);
	    return (role == REMOTE_DB || 
		    role == REMOTE_FUNC_FILES ||
		    role == REMOTE_DATA_PKGS);
	}
    } else if (item->flag == BTN_EXEC || item->flag == BTN_ADD) {
	return (role == FUNC_FILES);
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
	} else if (role == REMOTE_ADDONS) {
	    item->func = G_CALLBACK(show_addon_info);
	}
    } else if (item->flag == BTN_INDX) {
	/* index: databases only */
	if (role == NATIVE_DB) {
	    item->func = G_CALLBACK(open_db_index);
	} else if (role == REMOTE_DB) {
	    item->func = G_CALLBACK(open_remote_db_index);
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
    }

    return (item->func != NULL);
}

static void make_files_toolbar (windata_t *vwin)
{
    GtkWidget *hbox, *button;
    GretlToolItem *item;
    int i;

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gretl_toolbar_new();

    for (i=0; i<n_files_items; i++) {
	item = &files_items[i];
	if (winlist_item(item)) {
	    vwin_toolbar_insert_winlist(vwin);
	} else if (files_item_get_callback(item, vwin->role)) {
	    button = gretl_toolbar_insert(vwin->mbar, item, item->func, vwin, -1);
	    if (item->flag == BTN_ADD) {
		g_object_set_data(G_OBJECT(vwin->mbar), "add-button", button);
		gtk_widget_set_sensitive(button, FALSE);
	    }
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    vwin_add_finder(vwin);
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

/* handle drag of pointer from remote database window */

static void  
db_window_handle_drag  (GtkWidget *widget,
			GdkDragContext *context,
			gint x,
			gint y,
			GtkSelectionData *data,
			guint info,
			guint time,
			gpointer p)
{
#ifdef MAC_NATIVE
    if (info == GRETL_REMOTE_DB_PTR && data != NULL) {
	const guchar *seldata = gtk_selection_data_get_data(data);

	install_file_from_server(NULL, *(void **) seldata);
    }
#else
    if (info == GRETL_REMOTE_DB_PTR && data != NULL) {
	GdkAtom type = gtk_selection_data_get_data_type(data);
	
	if (type == GDK_SELECTION_TYPE_INTEGER) {
	    const guchar *seldata = gtk_selection_data_get_data(data);

	    install_file_from_server(NULL, *(void **) seldata);
	}
    }
#endif
}

/* handle drag of pointer from remote function package window */

static void  
pkg_window_handle_drag  (GtkWidget *widget,
			 GdkDragContext *context,
			 gint x,
			 gint y,
			 GtkSelectionData *data,
			 guint info,
			 guint time,
			 gpointer p)
{
#ifdef MAC_NATIVE
    if (info == GRETL_REMOTE_FNPKG_PTR && data != NULL) {
	const guchar *seldata = gtk_selection_data_get_data(data);

	install_file_from_server(NULL, *(void **) seldata);
    }
#else
    if (info == GRETL_REMOTE_FNPKG_PTR && data != NULL) {
	GdkAtom type = gtk_selection_data_get_data_type(data);

	if (type == GDK_SELECTION_TYPE_INTEGER) {
	    const guchar *seldata = gtk_selection_data_get_data(data);

	    install_file_from_server(NULL, *(void **) seldata);
	}
    }
#endif
}

static void set_up_viewer_drag_target (windata_t *vwin)
{
    GCallback callback;
    int i;

    if (vwin->role == NATIVE_DB) {
	i = GRETL_REMOTE_DB_PTR;
	callback = G_CALLBACK(db_window_handle_drag);
    } else if (vwin->role == FUNC_FILES) {
	i = GRETL_REMOTE_FNPKG_PTR;
	callback = G_CALLBACK(pkg_window_handle_drag);
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

void display_files (int role, gpointer data)
{
    GtkWidget *filebox;
    windata_t *vwin;
    struct fpkg_response fresp;
    gchar *title = NULL;
    int err = 0;

    vwin = get_browser_for_role(role);
    if (vwin != NULL) {
	gtk_window_present(GTK_WINDOW(vwin->main));
	return;
    }

    if (role == FUNC_FILES || role == NATIVE_DB) {
	title = files_title(role);
    } else if (role == PS_FILES) {
	title = g_strdup(_("gretl: practice files"));
    } else if (role == TEXTBOOK_DATA) {
	title = g_strdup(_("gretl: data files"));
    } else if (role == REMOTE_DB) {
	title = g_strdup(_("gretl: databases on server"));
    } else if (role == REMOTE_FUNC_FILES) {
	title = g_strdup(_("gretl: function packages on server"));
    } else if (role == REMOTE_DATA_PKGS) {
	title = g_strdup(_("gretl: data packages on server"));
    } else if (role == REMOTE_ADDONS) {
	title = g_strdup(_("gretl: addons"));
    }

    vwin = gretl_browser_new(role, title);
    g_free(title);

    fpkg_response_init(&fresp, data);

    if (role == REMOTE_DB) {
	gtk_window_set_default_size(GTK_WINDOW(vwin->main), 640, 480);
    }

    /* vertical box to hold file-listing widget and other elements */
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    make_files_toolbar(vwin);

    if (role == TEXTBOOK_DATA || role == PS_FILES) {
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

	build_datafiles_popup(vwin);
	while ((collection = pop_file_collection(role))) {
	    g_signal_connect(G_OBJECT(collection->listbox), "button-press-event",
			     G_CALLBACK(popup_menu_handler), 
			     vwin->popup);
	}
	reset_files_stack(role);
    } else if (role == FUNC_FILES || role == REMOTE_FUNC_FILES ||
	       role == REMOTE_ADDONS) {
	g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
			 G_CALLBACK(funcfiles_popup_handler), 
			 vwin);
	if (role == FUNC_FILES) {
	    connect_menu_add_signal(vwin);
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
	vwin->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(vwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->status, FALSE, FALSE, 0);
    } 

    /* put stuff into list box(es) */
    if (role == TEXTBOOK_DATA || role == PS_FILES) {
	err = populate_notebook_filelists(vwin, filebox, role);
    } else if (role == FUNC_FILES) {
	err = populate_filelist(vwin, &fresp);
	if (!err && fresp.col1_width > 15) {
	    /* widen the file box if need be */
	    gint w, h;

	    gtk_widget_get_size_request(filebox, &w, &h);
	    w += 100;
	    gtk_widget_set_size_request(filebox, w, h);
	}
    } else if (role == NATIVE_DB) {
	gint w, h, ndb = 0;

	err = populate_filelist(vwin, &ndb);
	if (!err && ndb > 12) {
	    gtk_widget_get_size_request(filebox, &w, &h);
	    h += 100;
	    gtk_widget_set_size_request(filebox, w, h);
	}	    
    } else {
	err = populate_filelist(vwin, NULL);
    }

    if (err) {
	gtk_widget_destroy(vwin->main);
    } else {
	gtk_widget_show_all(vwin->main); 
	gtk_widget_grab_focus(vwin->listbox);
	if (role == NATIVE_DB || role == FUNC_FILES) {
	    set_up_viewer_drag_target(vwin);
	} 
    }

    if (err) {
	if (role == FUNC_FILES && fresp.try_server == 1) {
	    /* no function packages on local machine */
	    display_files(REMOTE_FUNC_FILES, data);
	} else {
	    return;
	}
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
    if (!strcmp(s, "SFAddons"))
	return REMOTE_ADDONS;
    return 0;
}

/* make a browser window to display a set of files: textbook
   data files, practice scripts, databases...  
*/

void show_files (GtkAction *action, gpointer p)
{
    int code = display_files_code(gtk_action_get_name(action));

    display_files(code, p);
}

void show_native_dbs (void)
{
    display_files(NATIVE_DB, NULL);
}

static int get_func_info (const char *path, char **pdesc, 
			  char **pver)
{
    int err = get_function_file_header(path, pdesc, pver);

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static int fn_file_is_duplicate (const char *fname, 
				 const char *version,
				 GtkListStore *store,
				 int imax)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *fname_i;
    gchar *version_i;
    int i, n, ret = 0;

    model = GTK_TREE_MODEL(store);

    if (imax == 0 || !gtk_tree_model_get_iter_first(model, &iter)) {
	return 0;
    }

    n = strlen(fname) - 4;

    for (i=0; i<imax; i++) {
	gtk_tree_model_get(model, &iter, 
			   0, &fname_i, 
			   1, &version_i,
			   -1);
	if (strncmp(fname, fname_i, n) == 0 &&
	    strcmp(version, version_i) == 0) {
	    ret = 1;
	} 
	g_free(fname_i);
	g_free(version_i);
	if (ret || !gtk_tree_model_iter_next(model, &iter)) {
	    break;
	}
    }

    return ret;
}

static int is_functions_dir (const char *path)
{
    int n = strlen(path) - 9;

    return n > 0 && !strcmp(path + n, "functions");
}

static int ok_gfn_path (const char *fullname, 
			const char *shortname,
			const char *dirname,
			GtkListStore *store, 
			GtkTreeIter *iter,
			int imax, 
			int *maxlen)
{
    char *descrip = NULL, *version = NULL;
    int err, ok = 0;

    err = get_func_info(fullname, &descrip, &version);

    if (!err && !fn_file_is_duplicate(shortname, version, store, imax)) {
	if (iter != NULL) {
	    gchar *fname = g_strdup(shortname);
	    int n = strlen(fname) - 4;

	    fname[n] = '\0'; /* chop off ".gfn" */
	    if (n > *maxlen) {
		*maxlen = n;
	    }
	    gtk_list_store_append(store, iter);
	    gtk_list_store_set(store, iter, 
			       0, fname, 
			       1, version,
			       2, descrip, 
			       3, function_package_is_loaded(fullname), 
			       4, dirname, -1);
	    g_free(fname);
	}
	ok = 1;
    }

    free(descrip);
    free(version);

    return ok;
}

/* Read (or just count) the .gfn files in a given directory.
   The signal to count rather than read is that the @iter
   argument is NULL.
*/

static void
read_fn_files_in_dir (DIR *dir, const char *path, 
		      GtkListStore *store, GtkTreeIter *iter,
		      int *nfn, int *maxlen)
{
    struct dirent *dirent;
    char fullname[MAXLEN];
    int imax = *nfn;

    /* look first for the new-style setup: gfn file in its own
       subdir, e.g. functions/gig/gig.gfn 
    */

    if (is_functions_dir(path)) {
	while ((dirent = readdir(dir)) != NULL) {
	    const char *basename = dirent->d_name;
	
	    if (!strcmp(basename, ".") ||
		!strcmp(basename, "..")) {
		continue;
	    }

	    build_path(fullname, path, basename, NULL);
	    if (gretl_isdir(fullname)) {
		gchar *realbase, *realpath;
		FILE *fp;

		strcat(fullname, SLASHSTR);
		strcat(fullname, basename);
		strcat(fullname, ".gfn");
		fp = gretl_fopen(fullname, "r");
		if (fp != NULL) {
		    fclose(fp);
		    realbase = g_strdup_printf("%s.gfn", basename);
		    realpath = g_strdup_printf("%s%c%s", path, SLASH, basename);
		    *nfn += ok_gfn_path(fullname, realbase, realpath,
					store, iter, imax, maxlen);
		    g_free(realbase);
		    g_free(realpath);
		} else {
		    gretl_error_clear();
		}	
	    }
	}
    
	imax = *nfn;
	rewinddir(dir);
    }

    /* then look for "plain" .gfn files */

    while ((dirent = readdir(dir)) != NULL) {
	const char *basename = dirent->d_name;

	if (!strcmp(basename, ".") ||
	    !strcmp(basename, "..")) {
	    continue;
	}

	if (has_suffix(basename, ".gfn")) {
	    build_path(fullname, path, basename, NULL);
	    *nfn += ok_gfn_path(fullname, basename, path,
				store, iter, imax, maxlen);
	}
    }
}

static int dirname_done (char **dnames, int ndirs, char *dirname)
{
    int i;

    for (i=0; i<ndirs; i++) {
	if (!strcmp(dnames[i], dirname)) {
	    return 1;
	}
    }

    return 0;
}

gint populate_func_list (windata_t *vwin, struct fpkg_response *fresp)
{
    GtkListStore *store;
    GtkTreeIter iter;
    char **dnames = NULL;
    DIR *dir;
    int i, n_dirs = 0;
    int nfn, maxlen = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    nfn = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(vwin->listbox), "nfn"));

    if (nfn > 0) {
	/* we're repopulating an existing list */
	GtkTreeModel *model = GTK_TREE_MODEL(store);
	gchar *dirname;

	gtk_tree_model_get_iter_first(model, &iter);
	gtk_tree_model_get(model, &iter, 4, &dirname, -1);
	strings_array_add(&dnames, &n_dirs, dirname);
	g_free(dirname);

	while (gtk_tree_model_iter_next(model, &iter)) {
	    gtk_tree_model_get(model, &iter, 4, &dirname, -1);
	    if (!dirname_done(dnames, n_dirs, dirname)) {
		strings_array_add(&dnames, &n_dirs, dirname);
	    }
	    g_free(dirname);
	}
    }

    nfn = 0;
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (n_dirs > 0) {
	for (i=0; i<n_dirs; i++) {
	    dir = gretl_opendir(dnames[i]);
	    if (dir != NULL) {
		read_fn_files_in_dir(dir, dnames[i], store, &iter,
				     &nfn, &maxlen);
		closedir(dir);
	    }
	}
	strings_array_free(dnames, n_dirs);
    } else {
	dnames = get_plausible_search_dirs(FUNCS_SEARCH, &n_dirs);

	for (i=0; i<n_dirs; i++) {
	    dir = gretl_opendir(dnames[i]);
	    if (dir != NULL) {
		read_fn_files_in_dir(dir, dnames[i], store, &iter,
				     &nfn, &maxlen);
		closedir(dir);
	    }
	}
	strings_array_free(dnames, n_dirs);
    }	    

    if (nfn == 0) {
	if (fresp->try_server == 0) {
	    int resp;

	    resp = yes_no_dialog(_("gretl: function packages"),
				 _("No gretl function packages were found on this computer.\n"
				   "Do you want to take a look on the gretl server?"),
				 0);
	    if (resp == GRETL_YES) {
		fresp->try_server = 1;
	    }
	} else {
	    warnbox(_("No gretl function packages were found on this computer."));
	}
	return 1;
    } else {
	g_object_set_data(G_OBJECT(vwin->listbox), "nfn", GINT_TO_POINTER(nfn));
	presort_treelist(vwin);
    }

    return 0;
}

static void revise_loaded_status (GtkWidget *listbox)
{
    GtkTreeModel *model;
    GtkListStore *store;
    GtkTreeIter iter;
    gchar *fullname;
    gchar *filename;
    gchar *filepath;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(listbox));
    store = GTK_LIST_STORE(model);

    if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter)) {
	return;
    }

    while (1) {
	gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 
			   0, &filename, 4, &filepath, -1);
	fullname = g_strdup_printf("%s%c%s.gfn", filepath, SLASH, 
				   filename);
	gtk_list_store_set(store, &iter, 
			   3, function_package_is_loaded(fullname), 
			   -1);
	g_free(filename);
	g_free(filepath);
	g_free(fullname);
	if (!gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter)) {
	    break;
	}
    }
}

/* update function package status after run, edit calls */

void maybe_update_func_files_window (int role)
{
    windata_t *vwin = get_browser_for_role(FUNC_FILES);

    if (vwin != NULL && vwin->listbox != NULL) {
	if (role == EDIT_FN_PKG) {
	    populate_func_list(vwin, NULL);
	} else {
	    revise_loaded_status(vwin->listbox);
	}
    }
}

gint populate_filelist (windata_t *vwin, gpointer p)
{
    if (vwin->role == NATIVE_DB) {
	return populate_dbfilelist(vwin, p);
    } else if (vwin->role == REMOTE_DB) {
	return populate_remote_db_list(vwin);
    } else if (vwin->role == REMOTE_FUNC_FILES) {
	return populate_remote_func_list(vwin);
    } else if (vwin->role == REMOTE_DATA_PKGS) {
	return populate_remote_data_pkg_list(vwin);
    } else if (vwin->role == FUNC_FILES) {
	return populate_func_list(vwin, p);
    } else if (vwin->role == REMOTE_ADDONS) {
	return populate_remote_addons_list(vwin);
    } else {
	return read_file_descriptions(vwin, p);
    }
}

static GtkWidget *files_vbox (windata_t *vwin) 
{
    const char *data_titles[] = {
	N_("File"), 
	N_("Summary")
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
    const char *func_titles[] = {
	N_("Package"), 
	N_("Version"),
	N_("Summary"), 
	N_("Loaded?")
    };
    const char *remote_func_titles[] = {
	N_("Package"), 
	N_("Version"),
	N_("Author"),
	N_("Summary"), 
	N_("Local status")
    };
    const char *addons_titles[] = {
	N_("Package"), 
	N_("Version"),
	N_("Date"),
	N_("Local status")
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
    GType func_types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_BOOLEAN,
	G_TYPE_STRING   /* hidden string: directory */
    };
    GType remote_func_types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_BOOLEAN  /* hidden boolean: zipfile? */
    };
    GType addons_types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING  /* hidden string: description */
    };
    const char **titles = data_titles;
    GType *types = types_2;
    int full_width = 500, file_height = 260;
    int hidden_col = 0;
    int use_tree = 0;
    GtkWidget *vbox;
    int cols = 2;

    switch (vwin->role) {
    case NATIVE_DB:
	titles = db_titles;
	cols = 3;
	hidden_col = TRUE;
	break;
    case REMOTE_DB:
	titles = remote_db_titles;
	cols = 3;
	full_width = 580;
	use_tree = 1;
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
	file_height = 300;
	break;
    case FUNC_FILES:
	titles = func_titles;
	types = func_types;
	cols = G_N_ELEMENTS(func_types);
	hidden_col = TRUE;
	full_width = 560;
	file_height = 320;
	break;
    case REMOTE_FUNC_FILES:
	titles = remote_func_titles;
	types = remote_func_types;
	cols = G_N_ELEMENTS(remote_func_types);
	hidden_col = TRUE;
	full_width = 580;
	file_height = 340;
	break;
    case REMOTE_ADDONS:
	titles = addons_titles;
	types = addons_types;
	cols = G_N_ELEMENTS(addons_types);
	hidden_col = TRUE;
	full_width = 400;
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
    vwin_add_list_box(vwin, GTK_BOX(vbox), cols, hidden_col, 
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
   collections (Ramanathan, Wooldridge, etc.) or practice scripts.
   The function creates the pages but does not yet fill them out.
*/

static GtkWidget *files_notebook (windata_t *vwin, int role)
{
    file_collection *collection;
    GtkWidget *notebook;
    GtkWidget *page;
    GtkWidget *label;
    int err = 0;

    if (role != TEXTBOOK_DATA && role != PS_FILES) {
	/* we shouldn't be here! */
	return NULL;
    }
    
    /* assemble the info we'll need */
    err = build_file_collections();
    if (err) {
	gui_errmsg(err);
	return NULL;
    }

    notebook = gtk_notebook_new();
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook), TRUE);

    while ((collection = pop_file_collection(role))) {
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
    }

    reset_files_stack(role);

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
    const char *title;
    int found = 0;
    int pgnum = 0;

    if (role == TEXTBOOK_DATA) {
	title = get_datapage();
    } else {
	title = get_scriptpage();
    }

    reset_files_stack(role);

    while ((collection = pop_file_collection(role))) {
	vwin->listbox = collection->listbox;
	populate_filelist(vwin, collection);
	if (*title != '\0' && !strcmp(collection->title, title)) {
	    selected = collection;
	    pgnum = found;
	}
	found++;
    }

    if (found == 0) {
	/* didn't find anything to show */
	return 1;
    }

    if (selected == NULL) {
	reset_files_stack(role);
	selected = pop_file_collection(role);
    }

    reset_files_stack(role);

    vwin->listbox = selected->listbox;
    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), pgnum);
    gtk_widget_grab_focus(vwin->listbox);

    return 0;
} 
