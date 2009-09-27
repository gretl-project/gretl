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

#include "gretl_xml.h"
#include "gretl_func.h"

#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>

#define N_BROWSER_TYPES 5
#define BROWSER_BUSY    1
#define BROWSER_OK      0

static GtkWidget *browsers[N_BROWSER_TYPES];

static GtkWidget *files_window (windata_t *fdata);
static GtkWidget *files_notebook (windata_t *fdata, int code);
static int populate_notebook_filelists (windata_t *fdata, 
					GtkWidget *notebook,
					int code);

typedef struct _file_collection file_collection;

struct _file_collection {
    char *path;
    char *descfile;
    char *title;
    int which;
    GtkWidget *page;
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
    RECOG_OK,
    RECOG_ERROR,
    RECOG_NOT
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
                          c == REMOTE_DATA_PKGS)

static void
read_fn_files_in_dir (int role, DIR *dir, const char *path, 
		      GtkListStore *store, GtkTreeIter *iter,
		      int *nfn, int *maxlen);

static void 
fpkg_response_init (struct fpkg_response *f, gpointer p)
{
    f->col1_width = 0;

    if (p == NULL || p == mdata) {
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

#if COLL_DEBUG
    fprintf(stderr, "full_path: got '%s' from '%s' + '%s'\n",
	    fpath, s1, s2);
#endif

    return fpath;
}

static int recognized_collection (file_collection *coll)
{
    const file_collection recognized_data[] = {
	{ "wooldridge", "jw_descriptions", "Wooldridge", COLL_DATA, NULL },
	{ "gujarati", "dg_descriptions", "Gujarati", COLL_DATA, NULL },
	{ "pwt56", "descriptions", "PWT 56", COLL_DATA, NULL },
	{ NULL, NULL, NULL, COLL_DATA, NULL }
    }; 
    const file_collection recognized_ps[] = {
	{ "pwt56", "ps_descriptions", "PWT 56", COLL_PS, NULL },
	{ NULL, NULL, NULL, COLL_PS, NULL }
    }; 
    int i;

    for (i=0; recognized_data[i].path != NULL; i++) {
	if (strstr(coll->path, recognized_data[i].path) &&
	    !strcmp(coll->descfile, recognized_data[i].descfile)) {
	    coll->title = malloc(strlen(recognized_data[i].title) + 1);
	    if (coll->title == NULL) {
		return RECOG_ERROR;
	    }
	    strcpy(coll->title, recognized_data[i].title);
	    coll->which = COLL_DATA;
	    return RECOG_OK;
	}
    }

    for (i=0; recognized_ps[i].path != NULL; i++) {
	if (strstr(coll->path, recognized_ps[i].path) &&
	    !strcmp(coll->descfile, recognized_ps[i].descfile)) {
	    coll->title = malloc(strlen(recognized_ps[i].title) + 1);
	    if (coll->title == NULL) {
		return RECOG_ERROR;
	    }
	    strcpy(coll->title, recognized_ps[i].title);
	    coll->which = COLL_PS;
	    return RECOG_OK;
	}
    }

    return RECOG_NOT;
} 

static int get_title_from_descfile (file_collection *coll)
{
    FILE *fp;
    char *test;
    char line[64], title[24];
    int err = 0;

    test = full_path(coll->path, coll->descfile);

    fp = gretl_fopen(test, "r");

    if (fp == NULL) {
	err = 1;
    } else if (fgets(line, sizeof line, fp) == NULL) {
	err = 1;
    } else {
	chopstr(line);

	if (sscanf(line, "# %23[^:]", title) != 1) {
	    err = 1;
	} else {
	    coll->title = g_strdup(title);
	    if (coll->title == NULL) {
		err = 1;
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
    if (coll->title != NULL) {
	free(coll->title);
    }
    free(coll);
}  

static file_collection *file_collection_new (const char *path,
					     const char *descfile,
					     int *err)
{
    file_collection *coll;
    int chk;

    *err = 0;

    coll = malloc(sizeof *coll);
    if (coll == NULL) {
	*err = 1;
	return NULL;
    }

    coll->title = NULL;
    
    coll->path = gretl_strdup(path);
    if (coll->path == NULL) {
	free(coll);
	*err = 1;
	return NULL;
    }

    coll->descfile = gretl_strdup(descfile);
    if (coll->descfile == NULL) {
	free(coll->path);
	free(coll);
	*err = 1;
	return NULL;
    }  

    if (strstr(coll->descfile, "ps_")) {
	coll->which = COLL_PS;
    } else {
	coll->which = COLL_DATA;
    }

    chk = recognized_collection(coll);

    if (chk == RECOG_ERROR) {
	free_file_collection(coll);
	*err = 1;
	coll = NULL;
    } else if (chk == RECOG_NOT) {
	chk = get_title_from_descfile(coll);
	if (chk) {
	    free_file_collection(coll);
	    coll = NULL;
	}
    }
    
    return coll;
}

static int compare_colls (const void *a, const void *b)
{
    const file_collection *ca = *(const file_collection **) a;
    const file_collection *cb = *(const file_collection **) b;
     
    return strcmp(ca->title, cb->title);
}

static void collection_stack_sort (file_collection **colls, int n)
{
    file_collection *tmp;
    int i;

    for (i=0; i<n; i++) {
	if (!strcmp(colls[i]->title, "Gretl")) {
	    if (i > 0) {
		tmp = colls[0];
		colls[0] = colls[i];
		colls[i] = tmp;
	    }
	    break;
	}
    }

    qsort(colls + 1, n - 1, sizeof *colls, compare_colls);
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

static file_collection *pop_data_collection (void)
{
    return collection_stack(NULL, STACK_POP_DATA);
}

static file_collection *pop_ps_collection (void)
{
    return collection_stack(NULL, STACK_POP_PS);
}

void destroy_file_collections (void)
{
    collection_stack(NULL, STACK_DESTROY);
}

static void reset_data_stack (void)
{
    collection_stack(NULL, STACK_RESET_DATA);
}

static void reset_ps_stack (void)
{
    collection_stack(NULL, STACK_RESET_PS);
}

static void sort_data_stack (void)
{
    collection_stack(NULL, STACK_SORT_DATA);
}

static void sort_ps_stack (void)
{
    collection_stack(NULL, STACK_SORT_PS);
}

static int get_file_collections_from_dir (const char *dname, DIR *dir)
{
    file_collection *coll;
    struct dirent *dirent;
    size_t n;
    int err = 0;

    while (!err && (dirent = readdir(dir))) { 
	if (strstr(dirent->d_name, "descriptions")) {
#if COLL_DEBUG
	    fprintf(stderr, "   %s: looking at '%s'\n", dname, dirent->d_name);
#endif
	    n = strlen(dirent->d_name);
	    if (!strcmp(dirent->d_name + n - 12, "descriptions")) {
		coll = file_collection_new(dname, dirent->d_name, &err);
		if (coll != NULL && !err) {
		    err = push_collection(coll);
		}
	    }
	}
    }

    return err;
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

static int seek_file_collections (int location)
{
    char *tmp = NULL;
    DIR *dir, *try;  
    struct dirent *dirent;
    char *subdir;
    int i = 0, err = 0;

    if (location == DATA_SEARCH) {
	tmp = g_strdup_printf("%sdata", gretl_home());
    } else if (location == SCRIPT_SEARCH) {
	tmp = g_strdup_printf("%sscripts", gretl_home());
    } else if (location == USER_SEARCH) {
	tmp = g_strdup(gretl_workdir());
	trim_slash(tmp);
    } else {
	return 1;
    }

 user_search_2:

    dir = opendir(tmp);
    if (dir == NULL) {
	return 1;
    }

#if COLL_DEBUG
    fprintf(stderr, "seeking file collections in '%s'\n", tmp);
#endif

    while (!err && (dirent = readdir(dir))) {
	if (!dont_go_there(dirent->d_name)) {
	    if (strcmp(dirent->d_name, ".")) {
		subdir = full_path(tmp, dirent->d_name);
	    } else {
		subdir = tmp;
	    }
	    try = opendir(subdir);
	    if (try != NULL) {
#if COLL_DEBUG
		fprintf(stderr, " trying in subdir '%s'\n", subdir);
#endif
		err = get_file_collections_from_dir(subdir, try);
#if COLL_DEBUG
		fprintf(stderr, " result: err = %d\n", err);
#endif
		closedir(try);
	    }
	}
    }

    closedir(dir);
    g_free(tmp);

    if (location == USER_SEARCH && i++ == 0) {
	tmp = gretl_default_workdir();
	if (tmp != NULL) {
	    goto user_search_2;
	}
    }

    return err;
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

static void print_data_collections (void)
{
    file_collection *coll;

    printf("\n*** Data collections:\n");
    while ((coll = pop_data_collection())) {
	print_collection(coll);
    }
    reset_data_stack();
}

static void print_script_collections (void)
{
    file_collection *coll;

    printf("\n*** Script collections:\n");

    while ((coll = pop_ps_collection())) {
	print_collection(coll);
    }
    reset_ps_stack();
}
#endif

static int build_file_collections (void)
{
    static int built;
    static int err;

    if (!built) {
	err = seek_file_collections(DATA_SEARCH);
	if (!err) {
	    err = seek_file_collections(SCRIPT_SEARCH);
	}
	if (!err) {
	    err = seek_file_collections(USER_SEARCH);
	}
	if (!err) {
	    sort_data_stack();
	    sort_ps_stack();
	}
	built = 1;
    }

#if COLL_DEBUG
    print_data_collections();
    print_script_collections();
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
	     !strcmp(p, ".bn7"))) {
	    *p = '\0';
	}
    }

    return s;
}

static int read_file_descriptions (windata_t *win, gpointer p)
{
    FILE *fp;
    GtkListStore *store;
    GtkTreeSelection *selection;
    GtkTreeIter iter;
    char line[MAXLEN];
    char *index;
    file_collection *coll = (file_collection *) p;

    index = full_path(coll->path, coll->descfile);

    fp = gretl_fopen(index, "r");
    if (fp == NULL) {
	return 1;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(win->listbox)));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    
    while (fgets(line, sizeof line, fp)) {
	char fname[24], descrip[80], data[64];

	if (*line == '#') continue;

	if (win->role == TEXTBOOK_DATA) {
	    if (sscanf(line, " \"%23[^\"]\",\"%79[^\"]\"", 
		       fname, descrip) == 2) {
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter, 
				   0, strip_extension(fname), 
				   1, descrip, -1);
	    }
	} else { /* script files */
	    if (sscanf(line, " \"%23[^\"]\",\"%79[^\"]\",\"%63[^\"]\"", 
		       fname, descrip, data) == 3) {
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter, 
				   0, strip_extension(fname), 
				   1, descrip, 
				   2, data, -1);
	    }
	}
    }

    fclose(fp);

    /* select the first row */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(win->listbox));
    gtk_tree_selection_select_iter(selection, &iter);
    
    return 0;
}

static void show_datafile_info (GtkWidget *w, gpointer data)
{
    char hdrname[MAXLEN];
    windata_t *vwin = (windata_t *) data;
    char *descrip;
    PRN *prn;
    gchar *fname;
    file_collection *coll;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &fname);
    coll = g_object_get_data(G_OBJECT(vwin->listbox), "coll");

    build_path(hdrname, coll->path, fname, ".gdt");
    g_free(fname);

    descrip = gretl_get_gdt_description(hdrname);

    if (descrip != NULL) {
	prn = gretl_print_new_with_buffer(descrip);
	view_buffer(prn, 80, 320, _("gretl: data header"), INFO, NULL);
    } else {
	errbox(_("Failed to retrieve description of data"));
	fprintf(stderr, I_("didn't get description from %s\n"), hdrname);
    }
}

void browser_open_data (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *datname;
    file_collection *coll;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &datname);
    coll = g_object_get_data(G_OBJECT(vwin->listbox), "coll");

    build_path(tryfile, coll->path, datname, ".gdt");
    g_free(datname);

    set_datapage(coll->title);

    verify_open_data(vwin, OPEN_DATA);
}

void browser_open_ps (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *fname;
    file_collection *coll;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &fname);
    coll = g_object_get_data(G_OBJECT(vwin->listbox), "coll");

    build_path(scriptfile, coll->path, fname, ".inp");
    g_free(fname);

    /* close the calling window */
    gtk_widget_destroy(GTK_WIDGET(vwin->main));

    set_scriptpage(coll->title);

    view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT);
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
    int err = load_function_package_from_file(fname);

    if (err) {
	gui_errmsg(err);
    } 

    return err;
}

static int gui_delete_fn_pkg (const char *fname, windata_t *vwin)
{
    char *msg = g_strdup_printf(_("Function package %s"), fname);
    const char *opts[] = {
	N_("Delete package file"),
	N_("Unload member functions")
    };
    int active[] = {1, 1};
    int resp, err = 0;

    resp = checks_only_dialog ("gretl", msg,
			       opts, 2, active, 0);
    g_free(msg);

    if (resp < 0 || (active[0] == 0 && active[1] == 0)) {
	/* canceled, or equivalent */
	return 0;
    }

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

void set_funcs_dir_callback (windata_t *vwin, char *path)
{
    DIR *dir = opendir(path);

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

	resp = radio_dialog("gretl", NULL, opts, 2, 0, 0);
	if (resp < 0) {
	    return;
	} else if (resp > 0) {
	    gtk_list_store_clear(store);
	    nfn = nfn0 = 0;
	}

	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	read_fn_files_in_dir(vwin->role, dir, path, store, &iter,
			     &nfn, &maxlen);
	if (nfn > nfn0) {
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
    fprintf(stderr, "browser_functions_handler: path='%s'\n", path);
#endif

    if (task == LOAD_FN_PKG) {
	err = gui_load_user_functions(path);
    } else if (task == DELETE_FN_PKG) {
	err = gui_delete_fn_pkg(path, vwin);
    } else if (task == VIEW_FN_PKG_INFO) {
	display_function_package_data(pkgname, path, VIEW_PKG_INFO);
    } else if (task == VIEW_FN_PKG_CODE) {
	display_function_package_data(pkgname, path, VIEW_PKG_CODE);
    } else if (task == EDIT_FN_PKG) {
	edit_function_package(path);
    } else if (task == CALL_FN_PKG) {
	call_function_package(path, vwin->main, &err);
    }

    g_free(pkgname);

    if (!err && task == CALL_FN_PKG) {
	maybe_update_func_files_window(task);
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

static void set_browser_status (windata_t *vwin, int status)
{
    int i = vwin->role - TEXTBOOK_DATA;

    if (status == BROWSER_BUSY) {
	browsers[i] = vwin->main;
    } else {
	browsers[i] = NULL;
    }
} 

static void free_browser (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    set_browser_status(vwin, BROWSER_OK);
    free_windata(NULL, data);
}

static gpointer get_browser_ptr (int role)
{
    int i = role - TEXTBOOK_DATA;

    return (gpointer) &(browsers[i]);
}

static gpointer get_browser (int role)
{
    int i = role - TEXTBOOK_DATA;

    return browsers[i];
}

static int browser_busy (guint code)
{
    int ret = 0;

    if (code >= TEXTBOOK_DATA && code <= REMOTE_DB) {
	GtkWidget *w = browsers[code - TEXTBOOK_DATA];

	if (w != NULL) {
	    gtk_window_present(GTK_WINDOW(w));
	    ret = 1;
	}
    }

    return ret;
}

windata_t *get_local_viewer (int remote_role)
{
    windata_t *vwin = NULL;
    GtkWidget *w = NULL;

    if (remote_role == REMOTE_DB) {
	w = get_browser(NATIVE_DB);
    } else if (remote_role == REMOTE_FUNC_FILES) {
	w = get_browser(FUNC_FILES);
    }

    if (w != NULL) {
	vwin = g_object_get_data(G_OBJECT(w), "vwin");
    }

    return vwin;
}

static void new_package_callback (GtkWidget *w , gpointer p)
{
    data_save_selection_wrapper(SAVE_FUNCTIONS);
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

static void build_funcfiles_popup (windata_t *vwin)
{
    if (vwin->popup == NULL) {
	vwin->popup = gtk_menu_new();

	if (vwin->role == FUNC_FILES) {
	    /* local function files: full menu */
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
	    add_popup_item(_("Delete"), vwin->popup, 
			   G_CALLBACK(browser_del_func), 
			   vwin);
	    add_popup_item(_("New"), vwin->popup, 
			   G_CALLBACK(new_package_callback), 
			   vwin);
	} else {
	    /* files on server: limited menu */
	    add_popup_item(_("Info"), vwin->popup, 
			   G_CALLBACK(pkg_info_from_server), 
			   vwin);
	    add_popup_item(_("Install"), vwin->popup, 
			   G_CALLBACK(install_file_from_server), 
			   vwin);
	}
    }
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
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 6)
    { N_("Select directory"), GTK_STOCK_OPEN, NULL, BTN_DIR },
#else
    { N_("Select directory"), GTK_STOCK_DIRECTORY, NULL, BTN_DIR },
#endif
    { N_("Edit"),           GTK_STOCK_EDIT,       G_CALLBACK(browser_edit_func), BTN_EDIT },
    { N_("Info"),           GTK_STOCK_INFO,       NULL,                          BTN_INFO },
    { N_("View code"),      GTK_STOCK_PROPERTIES, G_CALLBACK(show_function_code), BTN_CODE },
    { N_("List series"),    GTK_STOCK_INDEX,      NULL,                          BTN_INDX },
    { N_("Install"),        GTK_STOCK_SAVE,       G_CALLBACK(install_file_from_server), BTN_INST },
    { N_("Execute"),        GTK_STOCK_EXECUTE,    G_CALLBACK(browser_call_func), BTN_EXEC },
    { N_("Delete"),         GTK_STOCK_DELETE,     G_CALLBACK(browser_del_func),  BTN_DEL },
    { N_("Look on server"), GTK_STOCK_NETWORK,    NULL,                          BTN_WWW },
    { N_("Local machine"),  GTK_STOCK_HOME,       NULL,                          BTN_HOME },
    { N_("New"),            GTK_STOCK_NEW,        G_CALLBACK(new_package_callback), BTN_NEW },
    { N_("Close"),          GTK_STOCK_CLOSE,      G_CALLBACK(close_files_viewer), BTN_CLOSE }
};

static int n_files_items = G_N_ELEMENTS(files_items);

#define common_item(f) (f == BTN_FIND || f == BTN_CLOSE)

#define local_funcs_item(f) (f == BTN_EDIT || f == BTN_NEW || \
			     f == BTN_DEL || f == BTN_CODE)

static int files_item_get_callback (GretlToolItem *item, int role)
{
    if (common_item(item->flag)) {
	return 1;
    } else if (local_funcs_item(item->flag)) {
	return (role == FUNC_FILES);
    } else if (item->flag == BTN_INST) {
	return (role == REMOTE_DB || 
		role == REMOTE_FUNC_FILES ||
		role == REMOTE_DATA_PKGS);
    } else if (item->flag == BTN_EXEC) {
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
    GtkWidget *hbox;
    GretlToolItem *item;
    int i;

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gretl_toolbar_new();

    for (i=0; i<n_files_items; i++) {
	item = &files_items[i];
	if (files_item_get_callback(item, vwin->role)) {
	    gretl_toolbar_insert(vwin->mbar, item, item->func, vwin, -1);
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    vwin_add_finder(vwin);
    gtk_widget_show_all(hbox);
}

#if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 8
/* we'll have to fake g_get_host_name */
static const gchar *g_get_host_name (void)
{
    return _("local machine");
}
#endif

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
    /* handle drag of pointer from remote database window */
    if (info == GRETL_REMOTE_DB_PTR && data != NULL && 
	data->type == GDK_SELECTION_TYPE_INTEGER) {
	install_file_from_server(NULL, *(void **) data->data);
    }
}

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
    /* handle drag of pointer from remote function package window */
    if (info == GRETL_REMOTE_FNPKG_PTR && data != NULL && 
	data->type == GDK_SELECTION_TYPE_INTEGER) {
	install_file_from_server(NULL, *(void **) data->data);
    }
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

void display_files (int code, gpointer p)
{
    GtkWidget *filebox;
    windata_t *vwin;
    struct fpkg_response fresp;
    gchar *title = NULL;
    int err = 0;

    if (browser_busy(code)) {
	/* an appropriate window is already open */
	return;
    }

    if (code == FUNC_FILES || code == NATIVE_DB) {
	title = files_title(code);
    } else if (code == PS_FILES) {
	title = g_strdup(_("gretl: practice files"));
    } else if (code == TEXTBOOK_DATA) {
	title = g_strdup(_("gretl: data files"));
    } else if (code == REMOTE_DB) {
	title = g_strdup(_("gretl: databases on server"));
    } else if (code == REMOTE_FUNC_FILES) {
	title = g_strdup(_("gretl: function packages on server"));
    } else if (code == REMOTE_DATA_PKGS) {
	title = g_strdup(_("gretl: data packages on server"));
    }

    vwin = gretl_browser_new(code, title, 0);
    g_free(title);

    fpkg_response_init(&fresp, p);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(free_browser),
		     vwin);

    set_browser_status(vwin, BROWSER_BUSY);

    if (code == REMOTE_DB) {
	gtk_window_set_default_size(GTK_WINDOW(vwin->main), 640, 480);
    }

    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    make_files_toolbar(vwin);

    if (code == TEXTBOOK_DATA || code == PS_FILES) {
	/* we'll need more than one tab */
	filebox = files_notebook(vwin, code);
    } else {
	filebox = files_window(vwin);
    }

    gtk_box_pack_start(GTK_BOX(vwin->vbox), filebox, TRUE, TRUE, 0);

    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);

    if (code == TEXTBOOK_DATA) { 
	file_collection *coll;

	build_datafiles_popup(vwin);
	while ((coll = pop_data_collection())) {
	    g_signal_connect(G_OBJECT(coll->page), "button-press-event",
			     G_CALLBACK(popup_menu_handler), 
			     vwin->popup);
	}
	reset_data_stack();
    } else if (code == FUNC_FILES || code == REMOTE_FUNC_FILES) {
	build_funcfiles_popup(vwin);
	g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
			 G_CALLBACK(popup_menu_handler), 
			 vwin->popup);
    } else if (code == NATIVE_DB || code == REMOTE_DB) {
	build_db_popup(vwin);
	g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
			 G_CALLBACK(popup_menu_handler), 
			 vwin->popup);
    } else if (code == REMOTE_DATA_PKGS)  {
	build_data_pkg_popup(vwin);
	g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
			 G_CALLBACK(popup_menu_handler), 
			 vwin->popup);
    }	

    if (REMOTE_ACTION(code)) {
	GtkWidget *hbox;

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
	vwin->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(vwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->status, FALSE, FALSE, 0);
    } 

    /* put stuff into list box(es) */
    if (code == TEXTBOOK_DATA || code == PS_FILES) {
	err = populate_notebook_filelists(vwin, filebox, code);
    } else if (code == FUNC_FILES) {
	err = populate_filelist(vwin, &fresp);
	if (fresp.col1_width > 15) {
	    /* widen the file box if need be */
	    gint w, h;

	    gtk_widget_get_size_request(filebox, &w, &h);
	    w += 100;
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
	if (code == NATIVE_DB || code == FUNC_FILES) {
	    set_up_viewer_drag_target(vwin);
	} 
    }

    if (err) {
	if (code == FUNC_FILES && fresp.try_server == 1) {
	    /* no function packages on local machine */
	    display_files(REMOTE_FUNC_FILES, p);
	} else {
	    return;
	}
    }

    if (code != TEXTBOOK_DATA && code != PS_FILES) {
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

static int get_func_info (const char *fullname, char **pdesc, 
			  char **pver)
{
    int err = get_function_file_header(fullname, pdesc, pver);

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

    if (!gtk_tree_model_get_iter_first(model, &iter)) {
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

#define pkg_is_loaded(n) (get_function_package_by_filename(n) != NULL)

static void
read_fn_files_in_dir (int role, DIR *dir, const char *path, 
		      GtkListStore *store, GtkTreeIter *iter,
		      int *nfn, int *maxlen)
{
    struct dirent *dirent;
    char fullname[MAXLEN];
    gchar *fname;
    int n, imax = *nfn;
    int err;

    while ((dirent = readdir(dir)) != NULL) {
	if (!strcmp(dirent->d_name, ".") ||
	    !strcmp(dirent->d_name, "..")) {
	    continue;
	}
	fname = g_strdup(dirent->d_name);
	n = strlen(fname);
	if (!g_ascii_strcasecmp(fname + n - 4, ".gfn")) {
	    char *descrip = NULL, *version = NULL;

	    build_path(fullname, path, fname, NULL);
	    err = get_func_info(fullname, &descrip, &version);

	    if (!err) {
		if (!fn_file_is_duplicate(fname, version, store, imax)) {
		    gtk_list_store_append(store, iter);
		    n -= 4;
		    fname[n] = '\0';
		    if (n > *maxlen) {
			*maxlen = n;
		    }
		    gtk_list_store_set(store, iter, 
				       0, fname, 
				       1, version,
				       2, descrip, 
				       3, pkg_is_loaded(fullname), 
				       4, path, -1);
		    *nfn += 1;
		}
		free(descrip);
		free(version);
	    }
	}
	g_free(fname);
    }
}

static gint 
compare_pkgnames (GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b,
		  gpointer p)
{
    gchar *t1, *t2;
    gint ret;

    gtk_tree_model_get(model, a, 0, &t1, -1);
    gtk_tree_model_get(model, b, 0, &t2, -1);

    lower(t1);
    lower(t2);

    ret = strcmp(t1, t2);

    g_free(t1);
    g_free(t2);

    return ret;    
}

static void sort_pkglist (windata_t *vwin)
{
    GtkTreeModel *model;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));

    gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(model), 
					 0, GTK_SORT_ASCENDING);
    gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(model), 0,
				    compare_pkgnames, NULL, NULL);
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
    int ndirs = 0;
    DIR *dir;
    int i, nfn, maxlen = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    nfn = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(vwin->listbox), "nfn"));

    if (nfn > 0) {
	/* we're repopulating an existing list */
	GtkTreeModel *model = GTK_TREE_MODEL(store);
	gchar *dirname;

	gtk_tree_model_get_iter_first(model, &iter);
	gtk_tree_model_get(model, &iter, 4, &dirname, -1);
	strings_array_add(&dnames, &ndirs, dirname);
	g_free(dirname);

	while (gtk_tree_model_iter_next(model, &iter)) {
	    gtk_tree_model_get(model, &iter, 4, &dirname, -1);
	    if (!dirname_done(dnames, ndirs, dirname)) {
		strings_array_add(&dnames, &ndirs, dirname);
	    }
	    g_free(dirname);
	}
    }

    nfn = 0;
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (ndirs > 0) {
	for (i=0; i<ndirs; i++) {
	    dir = opendir(dnames[i]);
	    if (dir != NULL) {
		read_fn_files_in_dir(vwin->role, dir, dnames[i], store, &iter,
				     &nfn, &maxlen);
		closedir(dir);
	    }
	}

	free_strings_array(dnames, ndirs);
    } else {
	char fndir[FILENAME_MAX];

	for (i=0; i<5; i++) {
	    *fndir = '\0';

	    if (i == 0) {
		/* pick up any function files in system dir */
		build_path(fndir, gretl_home(), "functions", NULL);
	    } else if (i == 1) {
		/* plus any function files in the user's working dir... */
		strcpy(fndir, gretl_workdir());
	    } else if (i == 2) {
		/* and in any "functions" subdir thereof */
		build_path(fndir, gretl_workdir(), "functions", NULL);
	    } else if (i == 3) {
		/* plus any in the user's dotdir */
		build_path(fndir, gretl_dotdir(), "functions", NULL);
	    } else if (i == 4) {
		/* plus any in the default working dir, if not already searched */
		char *tmp = gretl_default_workdir();

		if (tmp != NULL) {
		    build_path(fndir, tmp, "functions", NULL);
		    g_free(tmp);
		} 
	    }

	    if (*fndir != '\0') {
		dir = opendir(fndir);
		if (dir != NULL) {
		    read_fn_files_in_dir(vwin->role, dir, fndir, store, &iter,
					 &nfn, &maxlen);
		    closedir(dir);
		}
	    }
	}
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
	sort_pkglist(vwin);
    }

    return 0;
}

static void revise_loaded_status (GtkWidget *lbox)
{
    char fullname[MAXLEN];
    GtkTreeModel *model;
    GtkListStore *store;
    GtkTreeIter iter;
    gchar *fname;
    gchar *fdir;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(lbox));
    store = GTK_LIST_STORE(model);

    if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter)) {
	return;
    }

    while (1) {
	gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 0, &fname, 
			   4, &fdir, -1);
	sprintf(fullname, "%s%c%s.gfn", fdir, SLASH, fname);
	g_free(fname);
	g_free(fdir);
	gtk_list_store_set(store, &iter, 
			   3, pkg_is_loaded(fullname), 
			   -1);
	if (!gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter)) {
	    break;
	}
    }
}

/* update function package status after run, edit calls */

void maybe_update_func_files_window (int action)
{
    GtkWidget *w = get_browser(FUNC_FILES);
    windata_t *vwin = NULL;
    GtkWidget *lbox = NULL;

    if (w != NULL) {
	vwin = g_object_get_data(G_OBJECT(w), "vwin");
    }

    if (vwin != NULL) {
	lbox = vwin->listbox;
    }

    if (lbox != NULL) {
	if (action == EDIT_FN_PKG) {
	    populate_func_list(vwin, NULL);
	} else {
	    revise_loaded_status(lbox);
	}
    }
}

gint populate_filelist (windata_t *vwin, gpointer p)
{
    if (vwin->role == NATIVE_DB) {
	return populate_dbfilelist(vwin);
    } else if (vwin->role == REMOTE_DB) {
	return populate_remote_db_list(vwin);
    } else if (vwin->role == REMOTE_FUNC_FILES) {
	return populate_remote_func_list(vwin);
    } else if (vwin->role == REMOTE_DATA_PKGS) {
	return populate_remote_data_pkg_list(vwin);
    } else if (vwin->role == FUNC_FILES) {
	return populate_func_list(vwin, p);
    } else {
	return read_file_descriptions(vwin, p);
    }
}

static GtkWidget *files_window (windata_t *vwin) 
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
	N_("Summary"), 
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
    GType types_4[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING
    };
    GType types_5[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_BOOLEAN,
	G_TYPE_STRING
    };

    const char **titles = data_titles;
    GType *types = types_2;

    int full_width = 500, file_height = 260;
    int hidden_col = 0;
    int use_tree = 0;

    GtkWidget *box;
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
	cols = 5;
	types = types_5;
	hidden_col = TRUE;
	full_width = 560;
	file_height = 320;
	break;
    case REMOTE_FUNC_FILES:
	titles = remote_func_titles;
	cols = 4;
	types = types_4;
	full_width = 580;
	break;
    default:
	break;
    }

    if (cols == 3) {
	types = types_3;
    }

    full_width *= gui_scale;
    file_height *= gui_scale;

    box = gtk_vbox_new(FALSE, 0);
    gtk_widget_set_size_request(box, full_width, file_height);
    vwin_add_list_box(vwin, GTK_BOX(box), cols, hidden_col, 
		      types, titles, use_tree);
    gtk_widget_show(box);

    return box;
}

static void 
switch_file_page_callback (GtkNotebook *notebook, GtkNotebookPage *page,
			   guint page_num, windata_t *vwin)
{
    gpointer p = g_object_get_data(G_OBJECT(notebook), "browse_ptr");

    if (p != NULL) {
	GtkWidget *w = *(GtkWidget **) p;

	if (w != NULL) {
	    char wnum[4];

	    sprintf(wnum, "%d", (int) page_num);
	    vwin->listbox = g_object_get_data(G_OBJECT(notebook), wnum);
	}
    }
}

/* below: Construct a set of notebook pages for either data files
   (Ramanathan, Wooldridge, etc.) or practice scripts.  Creates the
   pages but does not yet fill them out.
*/

static GtkWidget *files_notebook (windata_t *vwin, int code)
{
    GtkWidget *notebook;
    GtkWidget *listpage;
    GtkWidget *label;
    char wnum[4];
    file_collection *coll;
    int j;

    if (code != TEXTBOOK_DATA && code != PS_FILES) {
	return NULL;
    }

    build_file_collections(); /* FIXME check for errors */

    notebook = gtk_notebook_new();
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook), TRUE);

    j = 0;

    while (1) {
	if (code == TEXTBOOK_DATA) {
	    coll = pop_data_collection();
	} else {
	    coll = pop_ps_collection();
	}

	if (coll == NULL) break;

	listpage = files_window(vwin);
	label = gtk_label_new(coll->title);
	gtk_widget_show(label);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), listpage, label);
	coll->page = vwin->listbox;
	sprintf(wnum, "%d", j);
	g_object_set_data(G_OBJECT(notebook), wnum, coll->page);
	g_object_set_data(G_OBJECT(coll->page), "coll", coll);
	g_signal_connect(G_OBJECT(coll->page), "key-press-event",
			 G_CALLBACK(enter_opens_file), vwin);
	j++;
    }

    if (code == TEXTBOOK_DATA) {
	reset_data_stack();
    } else {
	reset_ps_stack();
    }

    g_object_set_data(G_OBJECT(GTK_NOTEBOOK(notebook)), "browse_ptr",
		      get_browser_ptr(vwin->role));
    g_signal_connect(G_OBJECT(GTK_NOTEBOOK(notebook)), "switch-page",
		     G_CALLBACK(switch_file_page_callback),
		     vwin);

    if (gtk_notebook_get_n_pages(GTK_NOTEBOOK(notebook)) > 5) {
	gtk_notebook_popup_enable(GTK_NOTEBOOK(notebook));
    }

    gtk_widget_show(notebook);

    return notebook;
}

/* below: fill out a set of notebook pages (for data files
   or scripts files), entering the details into the page.
*/

static int populate_notebook_filelists (windata_t *vwin, 
					GtkWidget *notebook,
					int code)
{
    file_collection *coll;
    const char *title;
    int gotpref = 0;
    int gotcol = 0;
    int j;

    if (vwin == NULL) {
	return 1;
    }

    while (1) {
	if (code == TEXTBOOK_DATA) {
	    coll = pop_data_collection();
	} else {
	    coll = pop_ps_collection();
	}

	if (coll != NULL) {
	    gotcol = 1;
	} else {
	    break;
	}

	vwin->listbox = coll->page;
	populate_filelist(vwin, coll);
    }

    if (!gotcol) {
	return 1;
    }

    j = 0;

    if (code == TEXTBOOK_DATA) {
	reset_data_stack();
	title = get_datapage();
	if (*title != '\0') {
	    while ((coll = pop_data_collection())) {
		if (!strcmp(coll->title, title)) {
		    gotpref = 1;
		    break;
		}
		j++;
	    }
	    if (!gotpref) {
		reset_data_stack();
		coll = pop_data_collection();
		j = 0;
	    }
	} else {
	    coll = pop_data_collection();
	}
    } else {
	reset_ps_stack();
	title = get_scriptpage();
	if (*title != '\0') {
	    while ((coll = pop_ps_collection())) {
		if (!strcmp(coll->title, title)) {
		    gotpref = 1;
		    break;
		}
		j++;
	    }
	    if (!gotpref) {
		reset_ps_stack();
		coll = pop_ps_collection();
		j = 0;
	    }	    
	} else {
	    coll = pop_ps_collection();
	}
    }

    vwin->listbox = coll->page;

    if (code == TEXTBOOK_DATA) {
	reset_data_stack();
    } else {
	reset_ps_stack();
    }

    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), j);
    gtk_widget_grab_focus(vwin->listbox);

    return 0;
} 
